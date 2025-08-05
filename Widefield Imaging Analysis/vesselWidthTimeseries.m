function [vesselWidths, correlationCoefficients] = vesselWidthTimeseries(video,polyx,polyy,fitType,save_dir,save_filename)
%video: MxNxF single-channel video
%x: 5-point x coordinates of box 
%y: 5-point y coordinates of box 

meanImage = mean(video,3,'omitnan');

%calculate vessel diameters for each vessel box
imSize = size(video);
nf = imSize(3);
bf = 1;
op = 5;
nv = size(polyx,1);
vesselWidths = nan(nv, nf);
vesselProfiles = nan(nv, 1, 2); %[vessel, width(pixels), actualOrFit]
vesselEdges = nan(nv,2); %[vessel, left/right edge]
correlationCoefficients = nan(nv, nf);

waitfig = waitbar(0,'Calculating vessel widths');
for v = 1:nv
    waitbar((v-1)/nv,waitfig,'Calculating vessel widths')
    invert = false;

    %calculate vessel box position
    centerx = round(mean(polyx(v,1:4)));
    centery = round(mean(polyy(v,1:4)));

    %get vessel length and angle
    dx = mean(polyx(v,[2 3])) - mean(polyx(v,[1 4]));
    dy = mean(polyy(v,[2 3])) - mean(polyy(v,[1 4]));
    vesselLength = sqrt(dx^2 + dy^2);
    vesselAngle = atan2(dy,dx);

    %get vessel width
    dx = mean(polyx(v,[1 2])) - mean(polyx(v,[3 4]));
    dy = mean(polyy(v,[1 2])) - mean(polyy(v,[3 4]));
    vesselWidth = sqrt(dx^2 + dy^2);
    
    %get crop box around vessel
    minx = (imSize(2)/2) - (vesselLength/2);
    miny = (imSize(1)/2) - (vesselWidth);
    rect = round([minx miny vesselLength vesselWidth+(2*op)]);
    
    %store average vessel profile (actual and fit)
    vesselProfile = zeros(1, round(2*vesselWidth), 2);
    %store estimated vessel edge locations
    vesselEdge = zeros(1, 2);
    numavgs = 0;
    for f = 1:bf:nf
        waitbar(((v-1)/nv) + ((f-1)/(nf*nv)),waitfig,'Calculating vessel widths')

        if bf>1
            f_inds = f:(f+bf-1);
            f_inds(f_inds>nf) = [];
            I = nan([imSize(1:2) length(f_inds)]);
            for bfi = 1:length(f_inds)
                I(:,:,bfi) = video(:,:,f_inds(bfi));
            end
            I = mean(I,3,'omitnan');
        else
            I = video(:,:,f);
        end

        %center, rotate, and crop around vessel
        I = double(I);
        I = I - min(I,[],'all');
        I = I/max(I,[],'all');
        I = imtranslate(I, [((imSize(2)/2)-centerx), ((imSize(1)/2)-centery)]);  % center on vessel
        I = imrotate(I,rad2deg(vesselAngle),'crop');  % rotate so the vessel is horizontal
        I = imcrop(I,rect); %crop around vessel box

        %average vessel along length
        y = mean(I,2);
        y = y-min(y);
        y = y/max(y);
        nl = length(y);

        %if outside of vessel is brighter than the inside, invert brightness
        if f==1
            outside_brightness = mean(y([1:round(nl/4) round(3*nl/4):end]));
            inside_brightness = mean(y(round(nl/4):round(3*nl/4)));
            if outside_brightness > inside_brightness
                invert = true; %if the 1st frame gets inverted, invert all frames to stay consistent
            end
        end
        if invert
            y = 1-y; 
        end
        x = (1:nl)'; %create vessel profile vector (in pixels)

        % fit a curve to the vessel profile
        switch fitType
            case 'gaussian'
                gaussEqn = 'a*exp(-((x-b)/c)^2)';
                upper = [1.5 max(x) 0.3*length(x)];
                lower = [0.25 min(x) 0.01*length(x)];
                if f>1 && cc(1,2)>0.95 %if the previous loop found a good fit, use that as a starting point
                    startpoints = [fr.a fr.b fr.c];
                    fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',startpoints,'Upper',upper,'Lower',lower);
                else
                    startpoints = [1 mean(x) 0.1*length(x)];
                    fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',startpoints,'Upper',upper,'Lower',lower);
                end
                fr = fit(x,y,gaussEqn,fo);
                %evaluate fit
                yvec = fr.a*exp(-((x-fr.b)/fr.c).^2);
                %find edges by x-locations of left and right half-max of each gaussian term separately)
                left_edge = fr.b-fr.c*sqrt(-log(0.5));
                right_edge = fr.b+fr.c*sqrt(-log(0.5));

            case '2-term gaussian'
                if f>1 && cc(1,2)>0.95 %if the previous loop found a good fit, use that as a starting point
                    startpoints = [fr.a1 fr.b1 fr.c1 fr.a2 fr.b2 fr.c2];
                    fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',startpoints);
                    fr = fit(x,y,'gauss2',fo);
                else
                    fr = fit(x,y,'gauss2');
                end
                %evaluate fit
                yvec = fr.a1*exp(-((x-fr.b1)/fr.c1).^2) + fr.a2*exp(-((x-fr.b2)/fr.c2).^2);
                %find edges by x-locations of left and right half-max of each gaussian term separately)
                left_edge = min([fr.b1-fr.c1*sqrt(-log(0.5)), fr.b2-fr.c2*sqrt(-log(0.5))]);
                right_edge = max([fr.b1+fr.c1*sqrt(-log(0.5)), fr.b2+fr.c2*sqrt(-log(0.5))]);

            otherwise
                error('fit type not recognized')
        end
        
        %get correlation coefficient between profile and fit
        cc = corrcoef(y,yvec);
        correlationCoefficients(v,f) = cc(1,2);

        %store vessel profiles in average
        numavgs = numavgs+1;
        if nl>size(vesselProfile,2)
            vesselProfile(:,size(vesselProfile,2)+1:nl,:) = 0;
        end
        vesselProfile(1,1:nl,1) = (vesselProfile(1,1:nl,1)*(numavgs-1)/numavgs) + (y'/numavgs);
        vesselProfile(1,1:nl,2) = (vesselProfile(1,1:nl,2)*(numavgs-1)/numavgs) + (yvec'/numavgs);
        vesselWidths(v,f) = right_edge - left_edge;

        %store vessel edges in average
        vesselEdge = (vesselEdge*(numavgs-1)/numavgs) + ([left_edge right_edge]./numavgs);
    end
    vesselProfile_length = size(vesselProfile,2);
    if vesselProfile_length>size(vesselProfiles,2)
        vesselProfiles(:, size(vesselProfiles,2)+1:vesselProfile_length,:) = nan;
    end
    vesselProfiles(v, 1:vesselProfile_length, :) = vesselProfile;
    vesselEdges(v,:) = vesselEdge;
end

close(waitfig)

%fill missing values
if bf>1
    vesselWidths = fillmissing(vesselWidths,'spline',2,'MissingLocations',isnan(vesselWidths),'EndValues','nearest');
end

%get screen size for creating figure windows
screenSize = get(groot,'ScreenSize');
left = screenSize(3)*0.1;
bottom = screenSize(4)*0.1;
width = screenSize(3)*0.8;
height = 0.5*width;


%create a figure of all vessel gaussian fits
figure0 = figure('Position',[left bottom width height]);
npw = ceil(sqrt(nv));
nph = ceil(nv/npw);
for v = 1:nv
    subplot(nph,npw,v)
    plot(vesselProfiles(v,:,1),'k','LineWidth',2)
    hold on
    plot(vesselProfiles(v,:,2),'--b','LineWidth',2)
    plot([vesselEdges(v,1) vesselEdges(v,1)],[min(vesselProfiles(v,:,2)) max(vesselProfiles(v,:,2))],'g','LineWidth',2)
    plot([vesselEdges(v,2) vesselEdges(v,2)],[min(vesselProfiles(v,:,2)) max(vesselProfiles(v,:,2))],'g','LineWidth',2)
    title(['vessel ' num2str(v) ' (cc=' num2str(mean(correlationCoefficients(v,:))) ')'])
    xlabel('vessel cross section (pixels)')
    ylabel('brightness')
    if v==1
        legend({'mean vessel profile','mean gaussian fit','mean edges from fits'})
    end
end
saveas(figure0,fullfile(save_dir,['VesselFits_' save_filename '.svg']));

%create a figure of all drawn vessel boxes
figure1 = figure('Position',[left bottom width height]);
subplot(1,2,1)
axis1 = gca;
imshow((meanImage-0.05)*2,'Parent',axis1); %print mean image
set(figure1,'Position',[left bottom width height]);
set(axis1,'OuterPosition',[0 0 0.5 1]);
drawnow
hold(axis1,'on');
axis off
lw = 2;
colors = hsv(nv);
for v=1:nv %draw all ROIs
    plot(polyx(v,:),polyy(v,:),'Color',colors(v,:),'Parent',axis1,'LineWidth',lw);
end
pause(0.1)

%create a figure of all vessel diameter values
subplot(1,2,2)
axis2 = gca;
% axis2 = axes('Parent',figure2);
lw = 2;
colors = hsv(nv);

if nf==1
    %for single frame tif files, make scatterplot of
    %all vessels and diameters
    for v=1:nv %draw all ROIs
        scatter(v,vesselWidths(v),20,'filled','MarkerFaceColor',colors(v,:),'MarkerEdgeColor','none')
        hold on
    end
    xticks(1:nv)
    xlim([0.5 nv+0.5])
    ylabel('Vessel Width (pixels)')
    xlabel('Vessel Number')
    
else
    %for multi-frame tif files, draw line plots for all
    %vessel diameters
    for v=1:nv %draw all ROIs
        plot(vesselWidths(v,:),'Color',colors(v,:),'Parent',axis2,'LineWidth',lw);
        hold on
    end
    labels = cellstr(num2str((1:nv)'));
    legend(labels,'Parent',figure1,'Location','northeastoutside')
    ylabel('Vessel Width (pixels)')
    xlabel('Frame Number')
end