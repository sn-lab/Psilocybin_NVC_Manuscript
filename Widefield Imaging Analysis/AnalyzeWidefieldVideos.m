
%% set base directory
base_folder = 'D:\Widefield Data\New file organization\';


%% (STEP 1) get video directories
cond_names = {'before-control', 'after-control', 'before-psilocybin', 'after-psilocybin', 'before-doi', 'after-doi'};
wave_names = {'gcamp','hb','hbo','speckle','cmro'};

l = dir(base_folder);
for i = length(l):-1:1
    if l(i).name(1)=='.'
        l(i) = [];
    end
end

num_mice = length(l);
num_conds = length(cond_names);
num_wavelengths = length(wave_names);
foldernames = cell(num_mice,num_conds);
foundfiles = false(num_mice,num_conds);
for m = 1:num_mice
    for c = 1:num_conds
        cur_dir = fullfile(base_folder,l(m).name,cond_names{c},'Stim');
        msv4_dir = ls(fullfile(cur_dir,'MSV4*'));
        if ~isempty(msv4_dir)
            %check that vessel ROIs are present
            cur_dir = fullfile(cur_dir,msv4_dir);
            cur_name = ls(fullfile(cur_dir,'vesselWidths*.mat'));
            if isempty(cur_name)
                error(['cannot find vessel ROI .mat file for ' cur_dir])
            end
            foundfiles(m,c) = true;
            foldernames(m,c) = {cur_dir};
        end
    end
end

%delete files/folders in which no complete data was found
for m = num_mice:-1:1
    if all(foundfiles(m,:)==0)
        foundfiles(m,:) = [];
        foldernames(m,:) = [];
        l(m) = [];
    end
end
num_mice = length(l);

%print overall numbers
fprintf(['Total mice: ' num2str(num_mice) '\n'])
fprintf('Number of mice per condition: \n')
fprintf([num2str(sum(foundfiles)) '\n'])

% check for ROI numbers by vessel
roi_numbers = nan(num_mice,num_conds,2);
for m=1:num_mice
    for c = 1:num_conds
        if foundfiles(m,c)
            load(fullfile(foldernames{m,c},'vesselWidths_vessels_ref.mat'))
            roi_numbers(m,c,1) = sum(cellfun('length',vesselBoxes.type)==6); %num arteries
            roi_numbers(m,c,2) = sum(cellfun('length',vesselBoxes.type)==4); %num veins
        end
    end
end
max_rois = max(roi_numbers,[],'all');



%% (OPTIONAL) check that all foundfiles contain masks


for m=1:num_mice
    for c = 1:num_conds
        if foundfiles(m,c)
            if ~exist(fullfile(foldernames{m,c},'masks.mat'),"file")
                warning(['cannot find masks for mouse ' num2str(m) ', cond ' num2str(c)])
            end
        end
    end
end
fprintf('check for masks complete.\n')


%% (DONE) (run once) create masks for window, vessels, and gcamp labelling
% foldernames = {'D:\SN Lab\Psilocybin\Widefield Analysis\raw data\3108_Week1_before_stim\MSV4_3108_Week 1_before_stim\'};

start_mouse = 1;
start_cond = 1;
load_files = true;

for m = 1:num_mice
    if m>=start_mouse
        for c = 1:num_conds
            if c>=start_cond || m>start_mouse
                if foundfiles(m,c)
                    fprintf(['mouse ' num2str(m) ', cond ' num2str(c) '\n\n'])
            
                    %load files
                    if load_files
                        raw_ref = read_file(fullfile(foldernames{m,c},'raw_ref.tif'));
                        hbo = read_file(fullfile(foldernames{m,c},'hbo_mean.tif'));
                    end
                    
                    %%% set parameters
                    %window parameters:
                    prctile_range = [10 90]; %percentiles
                    medfilt_size = 75; %pixels width/height
                    range_threshold = 0.15; %binarization threshold
                    window_reduce_size = 0;
                    %vessel parameters:
                    binarization_sensitivity = 0.63;
                    open_size = 3; %size of small noise to remove
                    close_size = 5; %size of broken regions to connect
                    min_vessel_area = 200; %min size of vessel to keep
                    grow_size = 7; % # of pixels to grow
                    %gcamp parameters: 
                    gaussfilt_size = 25;
                    fraction_i530 = 0.25; %reduce the brightness of the i530 image
                    difference_threshold = 0.2; %fraction that gcamp has to be brighter than i530
                    gcamp_grow_size = 50; %number of pixels to grow around binarized region
                   
                    
                    figure()
                    %%% create mask of window region
                    %method: low range in hbo
                    i530 = double(raw_ref(:,:,3));
                    i530 = imgaussfilt(i530,3);
                    i530 = i530 - prctile(i530,1,'all');
                    i530 = i530/prctile(i530,99,'all');
                    subplot(3,4,1)
                    imshow(i530)
                    title('raw 530 nm')
                    
                    hbo_range = prctile(hbo,prctile_range(2),3) - prctile(hbo,prctile_range(1),3);
                    hbo_range = hbo_range - prctile(hbo_range,1,'all');
                    hbo_range = hbo_range/prctile(hbo_range,99,'all');
                    subplot(3,4,2)
                    imshow(hbo_range)
                    title('HbO 90-10 range')
                    
                    hbo_range = medfilt2(hbo_range,[75 75]);
                    subplot(3,4,3)
                    imshow(hbo_range)
                    title('HbO range filtered')
                    
                    windowMask = hbo_range>range_threshold;
                    se = strel('disk', window_reduce_size);
                    windowMask = ~imdilate(windowMask==0, se);
                    subplot(3,4,4)
                    imshow(windowMask)
                    title('window')
                    
                    
                    %%%create mask of vessels
                    % method: adaptive binarization
                    subplot(3,4,5)
                    imshow(i530)
                    title('raw 530 nm')
                    
                    binary_img = imbinarize(i530, 'adaptive', 'Sensitivity', binarization_sensitivity);
                    subplot(3,4,6)
                    imshow(binary_img)
                    title('adaptive binarization')
                    
                    se = strel('disk', open_size);
                    cleaned_img = imopen(~binary_img, se);
                    se = strel('disk', close_size);
                    cleaned_img = imclose(cleaned_img, se);
                    vesselMask = bwareaopen(cleaned_img, min_vessel_area);
                    subplot(3,4,7)
                    imshow(vesselMask & windowMask);
                    title('vessel & window')
                    
                    %grow dark region around vessels to create non-vessel mask
                    se = strel('disk', grow_size);
                    nonVesselMask = imdilate(vesselMask, se);
                    nonVesselMask = ~nonVesselMask;
                    subplot(3,4,8)
                    imshow(nonVesselMask & windowMask)
                    title('non-vessel & window')
                    
                    
                    %%% create mask of gcamp labelling region
                    %get gcamp image, only mask region
                    gcamp = double(raw_ref(:,:,2));
                    gcamp = imgaussfilt(gcamp,3);
                    gcamp = gcamp - prctile(gcamp,1,'all');
                    gcamp(~windowMask) = 0;
                    gcamp = gcamp/prctile(gcamp,99,'all');
                    subplot(3,4,9)
                    imshow(gcamp)
                    title('raw gcamp')
                    
                    %get i530 again, only mask region 
                    i530 = double(raw_ref(:,:,3));
                    i530 = imgaussfilt(i530,3);
                    i530 = i530 - prctile(i530,1,'all');
                    i530(~windowMask) = 0;
                    i530 = i530/prctile(i530,99,'all');

                    %gcamp region is relatively brighter than 530
                    gcamp = gcamp-(fraction_i530*i530);
                    subplot(3,4,10)
                    imshow(gcamp)
                    title('gcamp - i530')
                    
                    gcamp = imgaussfilt(double(gcamp),gaussfilt_size);
                    gcamp = gcamp>difference_threshold;
                    subplot(3,4,11)
                    imshow(gcamp)
                    title('filtered')
                    
                    se = strel('disk', gcamp_grow_size); % # of pixels to grow
                    gcamp = imdilate(gcamp, se);
                    gcampMask = gcamp;
                    subplot(3,4,12)
                    imshow(gcampMask & windowMask & nonVesselMask)
                    title('gcamp & window & non-vessel')
                    
                    answer = questdlg('Masks look OK?','manual evaluation of masks');
                    if strncmpi(answer,'y',1)
                        save(fullfile(foldernames{m,c},'masks.mat'),'windowMask','vesselMask','nonVesselMask','gcampMask')
                        clear raw_ref
                    else
                        error('set start mouse/condition and try again')
                    end
                end
            end
        end
    end
end



%% (OPTIONAL) view a specific mask

clear windowMask vesselMask nonVesselMask gcampMask
m = 1;
c = 1;
load(fullfile(foldernames{m,c},'masks.mat'))
I = read_file(fullfile(foldernames{m,c},'raw_ref.tif'));
figure()
subplot(2,2,1)
tmp = I(:,:,3);
tmp = tmp-prctile(tmp,1,'all'); tmp = tmp/prctile(tmp,75,'all');
imshow(tmp)
title('vessels')
subplot(2,2,2)
imshow(windowMask & nonVesselMask)
title('window')
subplot(2,2,3)
tmp = I(:,:,2);
tmp = tmp-prctile(tmp,1,'all'); tmp = tmp/prctile(tmp,95,'all');
imshow(tmp)
title('gcamp')
subplot(2,2,4)
imshow(gcampMask & windowMask & nonVesselMask)
title('gcamp & window & non-vessel')




%% (DONE) (run once) collapse videos into a single traces for each ROI

wb = waitbar(0,'');
for m = 1:num_mice
    for c = 1:num_conds
        waitbar((m-1)/num_mice + (c-1)/(num_conds*num_mice),wb,['loading mouse ' num2str(m), ', condition ' num2str(c)]);
        if foundfiles(m,c)
            if ~exist(fullfile(foldernames{m,c},'vessel_vectors.mat'),"file")
                %initialize space for vessel vectors
                vessel_vectors = nan(num_wavelengths,max_rois,3,200);
                %vessel_vectors: [wavelength(5), roi(6), roitype(1=Art,2=Vein,3=Tissue), frame(200)]

                for w = 1:num_wavelengths
                    %load file
                    I = read_file(fullfile(foldernames{m,c},[wave_names{w} '_mean.tif']));
                    [H,W,nf] = size(I);
    
                    %%%load vessel ROIs
                    load(fullfile(foldernames{m,c},'vesselWidths_vessels_ref.mat'))
                    for vtype = 1:2 %artery or vein
                        %get number of vessels of this type
                        if vtype==1 %artery
                            v_roi_inds = find(cellfun('length',vesselBoxes.type)==6);
                        else %vein
                            v_roi_inds = find(cellfun('length',vesselBoxes.type)==4);
                        end
                        num_v = length(v_roi_inds);

                        for v = 1:num_v
                            %create mask for this vessel
                            x = vesselBoxes.x(v_roi_inds(v),:);
                            y = vesselBoxes.y(v_roi_inds(v),:);
                            segmentMask = poly2mask(x,y,H,W);

                            % get mean value inside mask
                            for f = 1:nf
                                tmp = I(:,:,f);
                                tmp(isinf(tmp)) = nan;
                                vessel_vectors(w,v,vtype,f) = mean(tmp(segmentMask),'omitnan');
                            end
                        end
                    end
    
                    %%%load masks
                    load(fullfile(foldernames{m,c},'masks.mat'))
                    for f = 1:nf
                        tmp = I(:,:,f);
                        tmp(isinf(tmp)) = nan;
                        vessel_vectors(w,1,3,f) = mean(tmp(windowMask & nonVesselMask),'omitnan');
                        vessel_vectors(w,2,3,f) = mean(tmp(windowMask & nonVesselMask & gcampMask),'omitnan');
                    end
                end

                %save data for this mouse/condition
                save(fullfile(foldernames{m,c},'vessel_vectors.mat'),'vessel_vectors');
            end
        end
    end
end
close(wb)



%% (DONE) (run once) calculate vessel widths
start_mouse = 1;
start_cond = 1;
load_files = true;

wb = waitbar(0,'');
for m = 1:num_mice
    if m>=start_mouse
        for c = 1:num_conds
            if c>=start_cond || m>start_mouse

                waitbar((m-1)/num_mice + (c-1)/(num_conds*num_mice),wb,['loading mouse ' num2str(m), ', condition ' num2str(c)]);
                if foundfiles(m,c)
                    if ~exist(fullfile(foldernames{m,c},'vesselWidths_raw530_mean.mat'),"file")
                        if load_files
                            %load 530_mean video
                            video = read_file(fullfile(foldernames{m,c},'raw530_mean.tif'));
                
                            %load vesselBoxes
                            load(fullfile(foldernames{m,c},'vesselWidths_vessels_ref.mat'))
                        end
            
                        %calculate vessel widths across video
                        [vesselWidths, correlationCoefficients] = vesselWidthTimeseries(video,vesselBoxes.x,vesselBoxes.y,'gaussian',foldernames{m,c},'raw530_mean');
            
                        answer = questdlg('Fit look OK?','manual evaluation of profile fit');
                        if strncmpi(answer,'y',1)
                            save(fullfile(foldernames{m,c},'vesselWidths_raw530_mean.mat'),'vesselBoxes','vesselWidths','correlationCoefficients');
                            clear video
                        else
                            error('set start mouse/condition and try again')
                        end
                        save(fullfile(foldernames{m,c},'vesselWidths_raw530_mean.mat'),'vesselBoxes','vesselWidths','correlationCoefficients');
                    end
                end
            end
        end
    end
end
close(wb)


%% (OPTIONAL) spot-check vessel variables
m = 18;

for c = 1:6
    if foundfiles(m,c)
        load(fullfile(foldernames{m,c},'vessel_vectors.mat'),'vessel_vectors');
        
        % set timepoints of vectors
        fps = 10; %wasn't it resampled at 10 Hz?
        pretime = 3;
        t = ((1:200)/fps)-pretime;
        
        %plot vectors
        for mtype = 1:3
            figure('Position',[100 1100-(280*mtype) 700 250])
            for i = 1:5 %datatype
                for j = 1:6 %roi
                    subplot(5,6,j+(6*(i-1)))
                    y = squeeze(vessel_vectors(i,j,mtype,:));
                    x = t;
                    x(isnan(y)) = [];
                    y(isnan(y)) = [];
                    plot(x,y)
                    axis off
                end
            end
        end
        
        % plot widths
        figure('Position',[100 100 1200 150])
        load(fullfile(foldernames{m,c},'vesselWidths_raw530_mean.mat'));
        %vesselWidths
        t = t(1:size(vesselWidths,2));
        for v = 1:vesselBoxes.num
            subplot(1,12,v)
            switch vesselBoxes.type{v}
                case 'Artery'
                    plot(t,vesselWidths(v,:),'b')
                case 'Vein'
                    plot(t,vesselWidths(v,:),'r')
                otherwise
                    plot(t,vesselWidths(v,:),'k')
            end
            axis off
        end
        title(num2str(c))
    end
end


%% (OPTIONAL) inspect a specific vector
% vesselWidths(1,:)
% squeeze(vessel_vectors(5,1,3,:))'



%% (STEP 2) load previously processed mean ROI image vectors

%vesselwidths that didn't fit the profile well (exclude for now):
bad_widths = nan(23,5); %[mouse, vnum]
bad_widths(2,1:3) = [1 5 9];
bad_widths(3,1:2) = [7 10];
bad_widths(5,1) = 1;
bad_widths(6,1) = 8;
bad_widths(7,1:5) = [2 6 7 10 12];
bad_widths(8,1:2) = [2 10];
bad_widths(9,1) = 6;
bad_widths(10,1:5) = [3 4 7 9 12];
bad_widths(11,1:2) = [2 6];
bad_widths(12,1) = 10;
bad_widths(13,1:4) = [3 10 11 12];
bad_widths(14,1:3) = [4 7 11];
bad_widths(15,1:2) = [7 10];
bad_widths(16,1) = 8;
bad_widths(17,1) = 8;
bad_widths(18,1) = 6;
bad_widths(19,1:3) = [3 7 11]; %[noisy clipped noisy]
bad_widths(20,1) = 7; %clipped
bad_widths(21,1:3) = [7 9 11]; %[clipped noisy noisy]
bad_widths(22,1) = 7; %clipped
bad_widths(23,1) = 8; %clipped

all_vessel_vectors = nan(num_mice,num_conds,7,max_rois,3,200);
%all_vessel_vectors: [mouse, cond(6), wavelength(7), roi(6), roitype(1=Art,2=Vein,3=mask), frame(200)]
for m = 1:num_mice
    for c = 1:num_conds        
        if foundfiles(m,c)

            %load vessel vectors
            load(fullfile(foldernames{m,c},'vessel_vectors.mat'),'vessel_vectors');
            cur_max_rois = size(vessel_vectors,2);
            %vessel_vectors: [wavelength(5), roi(6), roitype(1=Art,2=Vein,3=mask), frame(200)]
            all_vessel_vectors(m,c,1:5,1:cur_max_rois,:,:) = vessel_vectors;

            %calculate HbT ("wavelength 6")
            all_vessel_vectors(m,c,6,:,:,:) = all_vessel_vectors(m,c,2,:,:,:) + all_vessel_vectors(m,c,3,:,:,:);

            %load vessel width data ("wavelength 7")
            load(fullfile(foldernames{m,c},'vesselWidths_raw530_mean.mat'))
            for vtype = 1:2 %artery or vein
                %get number of vessels of this type
                if vtype==1 %artery
                    v_roi_inds = find(cellfun('length',vesselBoxes.type)==6);
                else %vein
                    v_roi_inds = find(cellfun('length',vesselBoxes.type)==4);
                end
                num_v = length(v_roi_inds);
                nf = size(vesselWidths,2);
                
                for v = 1:num_v
                    %figure out what overall vessel # this is
                    if any(v_roi_inds(v)==bad_widths(m,:)) %bad vessel width - remove data
                        all_vessel_vectors(m,c,:,v,vtype,:) = nan;
                    else %good vessel fit - add width
                        % extract width vector for this vessel
                        tmp = vesselWidths(v_roi_inds(v),:);
                        all_vessel_vectors(m,c,7,v,vtype,1:nf) = tmp;
                    end
                end
            end
        end
    end
end


%% (STEP 3) process loaded data
% set timepoints of vectors
fps = 10;
pretime = 3;
t = ((1:200)/fps)-pretime;
baseline_frames = find(t<0);
stim_frames = find(t>0.5 & t<=6);
post_stim_frames = find(t>6.5 & t<=12);
stim_and_post_stim_frames = find(t>0.5 & t<=12);

[num_mice, num_conds, num_datatypes, max_vessels, num_vessel_types, num_frames] = size(all_vessel_vectors);
datanames = {'gcamp','hb','hbo','speckle','cmro','hbt','width'};
%vessel_vectors: [mouse, cond(6), datatype(7), roi(6), roitype(1=Art,2=Vein,3=mask), frame(200)]
%for roitypes 1 & 2, rois1-6 are vessels 1-6. for roitype 3, roi1 is parenchyma and roi2 is parenchyma-gcampregion

all_dvv_vectors = all_vessel_vectors;

% convert gcamp to % deltaX/X
tmp = all_dvv_vectors(:,:,1,:,:,:);
baseline_tmp = mean(tmp(:,:,:,:,:,baseline_frames),6,'omitnan');
baseline_tmp = repmat(baseline_tmp,[1 1 1 1 1 num_frames]);
tmp = tmp./baseline_tmp;
tmp = 100*(tmp-1);
all_dvv_vectors(:,:,1,:,:,:) = tmp;

% convert speckle to % deltaK/K
% tmp = all_dvv_vectors(:,:,4,:,:,:);
% baseline_tmp = mean(tmp(:,:,:,:,:,baseline_frames),6,'omitnan');
% baseline_tmp = repmat(baseline_tmp,[1 1 1 1 1 num_frames]);
% tmp = tmp./baseline_tmp;
% tmp = 100*(tmp-1);
% all_dvv_vectors(:,:,4,:,:,:) = tmp;

%calculate normalized ICT from cmro, hb, hbo
HbT0 = 100e-6;
HbR0 = 60e-6;
t_hb = all_dvv_vectors(:,:,2,:,:,:);
t_hbo = all_dvv_vectors(:,:,3,:,:,:);
t_cmro = all_dvv_vectors(:,:,5,:,:,:);
tmp = (1./(t_cmro+1)).*(1+(t_hb+t_hbo)./HbT0).*(1+t_hb/HbR0);
%convert ICT to % deltaS/S
baseline_tmp = mean(tmp(:,:,:,:,:,baseline_frames),6,'omitnan');
baseline_tmp = repmat(baseline_tmp,[1 1 1 1 1 num_frames]);
tmp = tmp./baseline_tmp;
tmp = 100*(tmp-1);
all_dvv_vectors(:,:,4,:,:,:) = tmp;

%convert CMRO to % deltaX/X (from very different baselines)
tmp = all_dvv_vectors(:,:,5,:,:,:);
baseline_tmp = mean(tmp(:,:,:,:,:,baseline_frames),6,'omitnan');
baseline_tmp = repmat(baseline_tmp,[1 1 1 1 1 num_frames]);
tmp = tmp./baseline_tmp;
tmp = 100*(tmp-1);
all_dvv_vectors(:,:,5,:,:,:) = tmp;

%convert Width to % deltaX/X (from very different baselines)
tmp = all_dvv_vectors(:,:,7,:,:,:);
baseline_tmp = mean(tmp(:,:,:,:,:,baseline_frames),6,'omitnan');
baseline_tmp = repmat(baseline_tmp,[1 1 1 1 1 num_frames]);
tmp = tmp./baseline_tmp;
tmp = 100*(tmp-1);
all_dvv_vectors(:,:,7,:,:,:) = tmp;

%create long format of vessel vectors
num_long_datatypes = 12;
all_long_vectors = nan(1000,num_long_datatypes,num_frames);
cur_row = 0;
for m = 1:num_mice
    for rt = 1:3
        for r = 1:max_rois
            for c = 1:num_conds
                tmp_row = nan(1,num_long_datatypes,num_frames);
                tmp_row(1,1,:) = all_dvv_vectors(m,c,1,r,rt,:); % 1. gcamp: % dF/F
                tmp_row(1,2,:) = all_dvv_vectors(m,c,2,r,rt,:); % 2. hb:
                tmp_row(1,3,:) = all_dvv_vectors(m,c,3,r,rt,:); % 3. hbo:
                tmp_row(1,4,:) = all_dvv_vectors(m,c,4,r,rt,:); % 4. cbf: % dS/S
                tmp_row(1,5,:) = all_dvv_vectors(m,c,5,r,rt,:); % 5. cmro:
                tmp_row(1,6,:) = all_dvv_vectors(m,c,6,r,rt,:); % 6. hbt:
                tmp_row(1,7,:) = all_dvv_vectors(m,c,7,r,rt,:); % 7. width: % dW/W
                tmp_row(1,8,:) = m; % 8. mouse ID
                tmp_row(1,9,:) = c; % 9. condition
                tmp_row(1,10,:) = rt; % 10. roitype: (1=artery, 2=vein, 3=parenchyma, 4=gcamp-parenchyma)
                tmp_row(1,11,:) = r; % 11. roi num
                tmp_row(1,12,:) = all_vessel_vectors(m,c,7,r,rt,:); % 12. width (pixels)
                if rt==3
                    if r==1
                        tmp_row(1,10,:) = 3; % 3=parenchyma
                    elseif r==2
                        tmp_row(1,10,:) = 4; % 4=gcamp-parenchyma)
                    end
                    tmp_row(1,11,:) = 1; %only 1 roi for this roitype
                end
                if any(~isnan(tmp_row(1,1:7,1))) %if there is any data for this mouse/cond/roi
                    cur_row = cur_row + 1;
                    all_long_vectors(cur_row,:,:) = tmp_row;
                end
            end
        end
    end
end

%average all vessels of the same type within each condition
num_vessels_per_condition = nan(num_mice,num_conds,2);
for m = 1:num_mice
    for c = 1:num_conds
        for vtype = 1:2
            num_vessels_per_condition(m,c,vtype) = sum(~isnan(all_dvv_vectors(m,c,1,:,vtype,1)));
        end
    end
end
all_dvv_vectors(:,:,:,1,1,:,:) = mean(all_dvv_vectors(:,:,:,:,1,:,:),4,'omitnan');
all_dvv_vectors(:,:,:,1,2,:,:) = mean(all_dvv_vectors(:,:,:,:,2,:,:),4,'omitnan');


%calculate mean response/mean recovery of all vectors
all_dvv_vals = mean(all_dvv_vectors(:,:,:,:,:,stim_frames),6,'omitnan');
all_dvv_vals(:,:,:,:,:,2) = mean(all_dvv_vectors(:,:,:,:,:,post_stim_frames),6,'omitnan');
all_long_vals = mean(all_long_vectors(:,:,stim_frames),3,'omitnan');
all_long_vals(:,:,2) = mean(all_long_vectors(:,:,post_stim_frames),3,'omitnan');

%notes
%all_dvv_vectors: [mouse, cond(6), datatype(7), mask(2, only for masks; 1 otherwise), roitype(1=Art,2=Vein,3=masks), frame(200)]
%all_dvv_vals: [mouse, cond(6), datatype(7), mask(2, only for masks; 1 otherwise), roitype(1=Art,2=Vein,3=masks), valtype(1=mean response,2=mean recovery)]
% 1. gcamp: % dF/F
% 2. hb:
% 3. hbo:
% 4. speed: % dS/S
% 5. cmro:
% 6. hbt:
% 7. width: % dW/W

save('C:\Users\misaa\Desktop\all_long_vectors.mat','all_long_vectors')

%% Set colors, figure sizes, etc

%SN vote winner red-to-blue palette
SN_colors = [0,   0,   0;...   %1.  black
             64,  64,  64;...  %2.  darkgrey
             128, 128, 128;... %3.  grey
             192, 192, 192;... %4.  lightgrey
             215, 48,  39;...  %5.  red
             252, 141, 89;...  %6.  pale-orange
             254, 224, 144;... %7.  light-pale-yellow
             224, 243, 248;... %8.  light-pale-blue
             145, 191, 219;... %9.  pale-blue
             69,  117, 180];   %10. blue
SN_colors = SN_colors/255; %RGB 0-1

brightwhite = [1 1 1];

colors = [SN_colors(3,:);...%pre control - grey
    SN_colors(1,:);...      %post control - black
    SN_colors(9,:);...      %pre psilo - pale-blue
    SN_colors(10,:);...     %post psilo - blue
    SN_colors(6,:);...      %pre mdl - pale-red
    SN_colors(5,:);...     %post mdl - red
    SN_colors(6,:);...      %pre psilo+mdl - pale-red
    SN_colors(5,:)];       %post psilo+mdl - red

alphas = [0.15 0.25 0.25 0.25 0.25 0.25 0.25 0.25];
linetypes = {'--','-','--','-','--','-','--','-'};
ylabels = {'ON neuron % \DeltaF/F','all vessel % \DeltaV/V','all vessel % \DeltaW/W'};
ylabels_baselines = {'raw brightness','velocity (\mm/s)','width (um)'};

vector_plot_size = [120 70];
box_plot_size = [120 85];
raster_plot_size = [130 250];
image_size = [85 85];
font_size = 4.8;


%% plot cartoon of camera/exposure timing
font_size = 4.8;
exposure_t = -0.1:0.001:1;
cam_t = zeros(size(exposure_t));
speckle_t = zeros(size(exposure_t));
gcamp_t = zeros(size(exposure_t));
i530_t = zeros(size(exposure_t));
i560_t = zeros(size(exposure_t));
trigger_dur = 0.05;
trigger_times = 0:0.1:1;
for i = 1:length(trigger_times)
    cam_t(exposure_t>=trigger_times(i) & exposure_t<(trigger_times(i)+trigger_dur)) = 1;
end
trigger_times = [0.1 0.5 0.9]-0.1;
for i = 1:length(trigger_times)
    speckle_t(exposure_t>=trigger_times(i) & exposure_t<(trigger_times(i)+trigger_dur)) = 1;
end
trigger_times = [0.2 0.6 1]-0.1;
for i = 1:length(trigger_times)
    gcamp_t(exposure_t>=trigger_times(i) & exposure_t<(trigger_times(i)+trigger_dur)) = 1;
end
trigger_times = [0.3 0.7]-0.1;
for i = 1:length(trigger_times)
    i530_t(exposure_t>=trigger_times(i) & exposure_t<(trigger_times(i)+trigger_dur)) = 1;
end
trigger_times = [0.4 0.8]-0.1;
for i = 1:length(trigger_times)
    i560_t(exposure_t>=trigger_times(i) & exposure_t<(trigger_times(i)+trigger_dur)) = 1;
end

figure('Position',[100 100 120 90])
axes('Position',[0.1 0.2 0.5 0.7])
plot(exposure_t,cam_t,'k')
hold on
plot(exposure_t,gcamp_t+1.5,'Color',[0 169 255]/255)
plot(exposure_t,i530_t+3,'Color',[94 255 0]/255)
plot(exposure_t,i560_t+4.5,'Color',[210 255 0]/255)
plot(exposure_t,speckle_t+6,'r')
text(1.1,0.5,'camera','Color','k','FontSize',font_size,'FontName','Arial')
text(1.1,2,'470 nm','Color',[0 169 255]/255,'FontSize',font_size,'FontName','Arial')
text(1.1,3.5,'530 nm','Color',[94 255 0]/255,'FontSize',font_size,'FontName','Arial')
text(1.1,5,'560 nm','Color',[210 255 0]/255,'FontSize',font_size,'FontName','Arial')
text(1.1,6.5,'785 nm','Color','r','FontSize',font_size,'FontName','Arial')
xlim([-0.05 1])
xlabel('time (s)')
ax = gca;
set(gca, 'YColor', 'none');
ylim([-0.5 7.5])
set(gca,'FontSize',font_size)
set(gca,'FontName','Arial')
box off
% print(gcf,'-vector','-dsvg',fullfile(base_folder,'exposure cartoon.svg'))


%% display example video frames and ROIs for each measurement
% data_folder = fullfile(base_folder,'raw data\MSV4_Less_habit_2R_Week 1_before_stim\');
data_folder = 'D:\SN Lab\Psilocybin\Widefield Analysis\raw data\MSV4_Less_habit_2R_Week 1_before_stim\';
image1names = {'speckle','gcamp','530','560'};
use_greyscale = true;

%create images of raw data channels surrounded by window mask
load(fullfile(data_folder,'masks.mat'))
load(fullfile(data_folder, 'vesselWidths_vessels_ref.mat'));

%get boundary line of window mask
boundaries = bwboundaries(windowMask);
boundary = boundaries{1};
boundary_x = boundary(:, 2);
boundary_y = boundary(:, 1);

% windowMask
% vesselMask
% nonVesselMask
% gcampMask
I = read_file(fullfile(data_folder,'raw_ref.tif'));
%speckle, gcamp, 530, 560
image_size = [68 68];
for i = 1:4
    tmp = I(:,:,i);
    tmp = tmp-prctile(tmp(:),2);
    tmp = tmp/prctile(tmp(:),98);
    if i>2
        tmp = tmp*3;
    end
    figure('Position',[100 100 image_size])
    ax = axes('Parent', gcf, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
    imshow(tmp)
    hold on
    plot(boundary_x,boundary_y,'--w','LineWidth',1)
    
    print(gcf,'-vector','-dsvg',fullfile(base_folder,['widefield image 1 ' image1names{i} '.svg']))

end

image_size = [85 85];


%%
%create images of normalized, stim-averaged files with gcamp/parenchyma/artery/vein masks
image2names = {'gcamp','ICT','530','hb','hbo','hbt'};
ylabels_colorbars = {'% \DeltaF/F','% \DeltaS/S','normalized 530 nm','normalized \DeltaHb','normalized \DeltaHbO','normalized \DeltaHbT'};
%show all in red/blue colormap of % deltaX/X

%gcamp
dt = 1;
I = read_file(fullfile(data_folder,'gcamp_mean.tif'));
IS = mean(I(:,:,stim_frames),3,'omitnan');
I0 = mean(I(:,:,baseline_frames),3,'omitnan');
half_range_prct = 15; % +/- range for colorbar
ID = (100/half_range_prct)*((IS./I0)-1);

figure('Position',[100 100 image_size])
ax = axes('Parent', gcf, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
if use_greyscale %with greyscale colormap
    ID(~windowMask) = 1;
    imshow(ID)
    hold on
    c = colorbar(gca,"eastoutside");
    c.Ticks = [0 1];
    c.TickLabels = {num2str(0) num2str(half_range_prct)};
else %with red/blur color map
    for i = 1:3
        tmp = interp1([-1 0 1],[SN_colors(10,i) 1 SN_colors(5,i)],ID,'linear');
        tmp(~windowMask) = 1;
        ID_rgb(:,:,i) = tmp;
    end
    imshow(ID_rgb)
    hold on
    map = interp1([-1 0 1],[0.2706 1 0.8431],-1:0.1:1,'linear');
    map = map';
    map(:,2) = interp1([-1 0 1],[0.4588 1 0.1882],-1:0.1:1,'linear');
    map(:,3) = interp1([-1 0 1],[0.7059 1 0.1529],-1:0.1:1,'linear');
    colormap(map)
    c = colorbar(gca,"eastoutside");
    c.Ticks = [0 1];
    c.TickLabels = {num2str(-half_range_prct) num2str(half_range_prct)};
end
c.Limits = [0 1];
print(gcf,'-vector','-dsvg',fullfile(base_folder,['widefield image 2 ' image2names{dt} '.svg']))

%draw ROI
if use_greyscale
    ID_rgb = repmat(ID,[1 1 3]);
end
for i = 1:3
    tmp = ID_rgb(:,:,i);
    tmp(windowMask & nonVesselMask & gcampMask) = SN_colors(7,i);
    ID_rgb(:,:,i) = tmp;
end

imshow(ID_rgb)
print(gcf,'-vector','-dsvg',fullfile(base_folder,['widefield image 2 ' image2names{dt} ' ROIs.svg']))

%%
% ICT
dt = 2;

HbT0 = 100e-6;
HbR0 = 60e-6;
cmro = read_file(fullfile(data_folder,'cmro_mean.tif'));
Hb = read_file(fullfile(data_folder,'hb_mean.tif'));
HbO = read_file(fullfile(data_folder,'hbo_mean.tif'));
I = (1./(cmro+1)).*(1+(Hb+HbO)./HbT0).*(1+Hb/HbR0);

IS = mean(I(:,:,stim_frames),3,'omitnan');
I0 = mean(I(:,:,baseline_frames),3,'omitnan');

half_range_prct = 20; % +/- range for colorbar


ID = (100/half_range_prct)*((IS./I0)-1);

figure('Position',[100 100 image_size])
ax = axes('Parent', gcf, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
if use_greyscale %with greyscale colormap
    ID(~windowMask) = 1;
    imshow(ID)
    hold on
    c = colorbar(gca,"eastoutside");
    c.Ticks = [0 1];
    c.TickLabels = {num2str(0) num2str(half_range_prct)};
else %with red/blur color map
    for i = 1:3
        tmp = interp1([-1 0 1],[SN_colors(10,i) 1 SN_colors(5,i)],ID,'linear');
        tmp(~windowMask) = 1;
        ID_rgb(:,:,i) = tmp;
    end
    imshow(ID_rgb)
    hold on
    map = interp1([-1 0 1],[0.2706 1 0.8431],-1:0.1:1,'linear');
    map = map';
    map(:,2) = interp1([-1 0 1],[0.4588 1 0.1882],-1:0.1:1,'linear');
    map(:,3) = interp1([-1 0 1],[0.7059 1 0.1529],-1:0.1:1,'linear');
    colormap(map)
    c = colorbar(gca,"eastoutside");
    c.Ticks = [0 1];
    c.TickLabels = {num2str(-half_range_prct) num2str(half_range_prct)};
end
c.Limits = [0 1];
print(gcf,'-vector','-dsvg',fullfile(base_folder,['widefield image 2 ' image2names{dt} '.svg']))

%draw ROI
if use_greyscale
    ID_rgb = repmat(ID,[1 1 3]);
end
for i = 1:3
    tmp = ID_rgb(:,:,i);
    tmp(windowMask & nonVesselMask) = SN_colors(7,i);
    ID_rgb(:,:,i) = tmp;
end
imshow(ID_rgb)
for i = 1:vesselBoxes.num
    switch(vesselBoxes.type{i})
        case 'Artery'
            patch(vesselBoxes.x(i,:),vesselBoxes.y(i,:),'k','EdgeColor',SN_colors(5,:),'FaceColor',SN_colors(5,:))
        case 'Vein'
            patch(vesselBoxes.x(i,:),vesselBoxes.y(i,:),'k','EdgeColor',SN_colors(10,:),'FaceColor',SN_colors(10,:))
    end
end
print(gcf,'-vector','-dsvg',fullfile(base_folder,['widefield image 2 ' image2names{dt} ' ROIs.svg']))


%%
% width
dt = 3;
I = read_file(fullfile(data_folder,'raw_ref.tif'));
I = I(:,:,3); %I530
ID = I;
ID = ID-prctile(ID(:),1);
ID = 3*ID/prctile(ID(:),99);

figure('Position',[100 100 image_size])
ax = axes('Parent', gcf, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
imshow(ID)
hold on

print(gcf,'-vector','-dsvg',fullfile(base_folder,['widefield image 2 ' image2names{dt} '.svg']))

%draw ROI
for i = 1:vesselBoxes.num
    switch(vesselBoxes.type{i})
        case 'Artery'
            plot(vesselBoxes.x(i,1:2),vesselBoxes.y(i,1:2),'Color',SN_colors(5,:),'LineWidth',1)
            plot(vesselBoxes.x(i,3:4),vesselBoxes.y(i,3:4),'Color',SN_colors(5,:),'LineWidth',1)
        case 'Vein'
            plot(vesselBoxes.x(i,1:2),vesselBoxes.y(i,1:2),'Color',SN_colors(10,:),'LineWidth',1)
            plot(vesselBoxes.x(i,3:4),vesselBoxes.y(i,3:4),'Color',SN_colors(10,:),'LineWidth',1)
    end
end
print(gcf,'-vector','-dsvg',fullfile(base_folder,['widefield image 2 ' image2names{dt} ' ROIs.svg']))



%%

% hb, hbo, hbt
hb = read_file(fullfile(data_folder,[image2names{4} '_mean.tif']));
hbo = read_file(fullfile(data_folder,[image2names{5} '_mean.tif']));
hbt = hb+hbo;

Imaxs = [10 25 25];
for dt = 4:6
    switch dt
        case 4
            I = hb;
        case 5
            I = hbo;
        case 6 
            I = hbt;
    end
    half_range_prct = 100; % +/- range for colorbar
    ID = mean(I(:,:,stim_frames),3,'omitnan')*1000000; %convert units
    % Imax = prctile(abs(ID(:)),99);
    Imax = Imaxs(dt-3);
    ID = (ID/Imax);
    ID(ID>1) = 1;
    ID(ID<-1) = -1;

    figure('Position',[100 100 image_size])
    ax = axes('Parent', gcf, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
    if use_greyscale %with greyscale colormap
        ID(~windowMask) = 1;
        if dt==4
            ID = ID/2 + 0.5; %shift so that negatives are visible
        end
        imshow(ID)
        c = colorbar(gca,"eastoutside");
        c.Ticks = [0 1];
        c.TickLabels = {num2str(0) num2str(half_range_prct)};
    else %with red/blur color map
        for i = 1:3
            tmp = interp1([-1 0 1],[SN_colors(10,i) 1 SN_colors(5,i)],ID,'linear');
            tmp(~windowMask) = 1;
            ID_rgb(:,:,i) = tmp;
        end
        imshow(ID_rgb)
    
        map = interp1([-1 0 1],[0.2706 1 0.8431],-1:0.1:1,'linear');
        map = map';
        map(:,2) = interp1([-1 0 1],[0.4588 1 0.1882],-1:0.1:1,'linear');
        map(:,3) = interp1([-1 0 1],[0.7059 1 0.1529],-1:0.1:1,'linear');
        colormap(map)
        c = colorbar(gca,"eastoutside");
        c.Ticks = [0 1];
        c.TickLabels = {num2str(-Imaxs(dt-3)) num2str(Imaxs(dt-3))};
    end

    c.Limits = [0 1];
    
    print(gcf,'-vector','-dsvg',fullfile(base_folder,['\widefield image 2 ' image2names{dt} '.svg']))

    %draw ROI
    if use_greyscale
        ID_rgb = repmat(ID,[1 1 3]);
    end
    for i = 1:3
        tmp = ID_rgb(:,:,i);
        tmp(windowMask & nonVesselMask) = SN_colors(7,i);
        ID_rgb(:,:,i) = tmp;
    end
    imshow(ID_rgb)
    for i = 1:vesselBoxes.num
        switch(vesselBoxes.type{i})
            case 'Artery'
                patch(vesselBoxes.x(i,:),vesselBoxes.y(i,:),'k','EdgeColor',SN_colors(5,:),'FaceColor',SN_colors(5,:))
            case 'Vein'
                patch(vesselBoxes.x(i,:),vesselBoxes.y(i,:),'k','EdgeColor',SN_colors(10,:),'FaceColor',SN_colors(10,:))
        end
    end
    print(gcf,'-vector','-dsvg',fullfile(base_folder,['widefield image 2 ' image2names{dt} ' ROIs.svg']))


    switch dt
        case 4
            tmpVenule = ID;
        case 6
            tmpArteriole = ID;
    end
end

% arteriole/venule map
AVmap = tmpArteriole;
AVmap(:,:,3) = -tmpVenule;
AVmap(:,:,2) = 0;
print(gcf,'-vector','-dsvg',fullfile(base_folder,'widefield avmap.svg'))


%% plot mean +/- SEM vectors for each condition

vector_plot_size = [120 75]; %smaller version

stimtime = [0 6];
xlims = [-3 15];
ylims = [-3 7; -10 10; -15 25; -8 17; -12 15; -15 25; -3 7]; %width as pixels
linetypes = {'--','-','--','-','--','-'};
num_mice_per_exp = nan(num_conds, num_datatypes, 4);
num_samples_per_exp = nan(num_conds, num_datatypes, 4);
% conds_to_plot = [2 1 4 3]; %control, psilocybin
conds_to_plot = [6 5]; %DOI
% 1. gcamp: % dF/F
% 2. hb:
% 3. hbo:
% 4. speed: % dS/S
% 5. cmro:
% 6. hbt:
% 7. width: % dW/W

%for BOLD modeling:
speed_means = nan(length(t),4,4); %[pts, roitype(4), cond(4)]

%plot all_dvv_vectors

rts = {4, [1 2 3],[1 2 3],[1 2 3],[],[1 2 3],[1 2]};
ylabels = {'arterioles','venules','parenchyma','parenchyma'};

for w = [1 2 3 4 6 7]
    for rt = rts{w}
        figure('Position',[100 100 vector_plot_size])

        plot([0 20],[0 0],'--k')
        hold on
        patch([0 6 6 0],[-100 -100 100 100],'k','EdgeColor','none','FaceAlpha',0.1)
        text(3,ylims(w,2),'stimulus','Color',[0.2 0.2 0.2],'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',font_size)

        ylineloc = ylims(w,1) + 0.05*(ylims(w,2)-ylims(w,1));
        ylineh = 0.02*(ylims(w,2)-ylims(w,1));
        if rt==rts{w}(1)
            if w==1 %plot response region
                plot([0.5 6],[ylineloc ylineloc],'k')
                plot([0.5 0.5],ylineloc+[ylineh -ylineh],'k')
                plot([6 6],ylineloc+[ylineh -ylineh],'k')
                text(3.25,ylineloc,'response','FontSize',font_size,'FontName','Arial','HorizontalAlignment','center','VerticalAlignment','bottom')
            else
                plot([0.5 2.5],[ylineloc ylineloc],'k')
                plot([2.5 12],[ylineloc ylineloc],'k')
                plot([0.5 0.5],ylineloc+[ylineh -ylineh],'k')
                plot([2.5 2.5],ylineloc+[ylineh -ylineh],'k')
                plot([12 12],ylineloc+[ylineh -ylineh],'k')
                text(1.5,ylineloc,'rise','FontSize',font_size,'FontName','Arial','HorizontalAlignment','center','VerticalAlignment','bottom')
                text(7.25,ylineloc,'recovery','FontSize',font_size,'FontName','Arial','HorizontalAlignment','center','VerticalAlignment','bottom')
            end
        end
            
        %plot vectors
        for c = conds_to_plot

            %find all rows in this condition, for this roitype, with non-nans in this datatype
            rows = find(all_long_vals(:,9,1)==c & all_long_vals(:,10,1)==rt & ~isnan(all_long_vals(:,w,1)));
            tmp = permute(all_long_vectors(rows,w,:),[1 3 2]);
            %tmp: [mouse, frame]
    
            num_mice_per_exp(c,w,rt) = length(unique(all_long_vals(rows,8,1)));
            num_samples_per_exp(c,w,rt) = length(rows);

            %get avg +/- sem of all mice for this condition
            tmp = fillmissing(tmp,"nearest",2,"MissingLocations",isnan(tmp));
            if w==2 || w==3 || w==6
                tmp = tmp*1000000;
            end
            meanvec = mean(tmp,'omitnan');
            semvec = std(tmp,'omitnan')/sqrt(num_samples_per_exp(c,w,rt));
            ul = meanvec+semvec;
            ll = meanvec-semvec;
            patch([t fliplr(t)],[ll fliplr(ul)],'k','EdgeColor','none','FaceColor',colors(c,:),'FaceAlpha',alphas(c))
            plot(t,meanvec,[linetypes{c} 'k'],'Color',colors(c,:),'LineWidth',1)

            if w==4
                speed_means(:,rt,c) = meanvec;
            end
        end
        box off
        xlim(xlims)
        ylim(ylims(w,:))
        xticks(0:5:15)
        box off
        xlabel('time (s)','FontSize',font_size)
        % ylabel(ylabels{w,rt},'FontSize',font_size)
        ylabel(ylabels{rt},'FontSize',font_size)
        set(gca,'FontSize',font_size)
        set(gca,'FontName','Arial')

        % legend({'a'},'Location','northeastoutside')
        print(gcf,'-vector','-dsvg',fullfile(base_folder,['widefield vector w' num2str(w) ' rt' num2str(rt) '.svg']))
    
        % close(gcf)
    end
end

%save artery and vein speeds means for BOLD simulation
% save('D:\SN Lab\Psilocybin\BOLD simulation\speed_avgs_widefield.mat','speed_means')


%% Box plots and stats
%neural activity box plots and stats

%all_long_vectors:
%[measurement, datatype, vector]
%datatypes:
%1. dF/F
%2. Hb
%3. HbO
%4. dS/S
%5. cmro
%6. hbt
%7. width
%8. mouseID
%9. condition
%10. roitype (1=artery, 2=vein, 3=parenchyma, 4=gcamp-parenchyma)
%11. roinum

%set time ranges for analysis of neuron/vessel response
neuron_response_inds = find(t>0.5 & t<=6);

do_stats = true;
do_plots = true;

box_plot_size = [120 60]; %smaller version
box_ylimits = [-6 5];
lme_treatment_ps = nan(6,3,2); %[dt, roitype, ]
cur_comparison = 0;
dtypenames = {'calcium','flow', 'hb','hbo', 'hbt','width'};
ylabels_boxplots = {'\Delta response (% \DeltaF/F)', '\Delta recovery (% \DeltaS/S)', '\Delta recovery (\DeltaHb)', '\Delta recovery (\DeltaHbO)', '\Delta recovery (\DeltaHbT)', '\Delta recovery (% \DeltaW/W)'};
curlabel = 0;
rt = 1; %neurons has only 1 roi type
dt = 1; %calcium

% conds_to_compare = [2 4];  %control pre/post, psilo pre/post
% xlabels_boxplots = {'ctrl','psil'};

conds_to_compare = [2 6];  %control pre/post, DOI pre/post
xlabels_boxplots = {'ctrl','DOI'};

%allocate space for delta values
pre_vals = nan(1,8);
post_vals = nan(1,8);
delta_vals = nan(1,8);
Data = [];
Mouse = [];
Treatment = [];
Timepoint = [];
Cell = [];

for c = conds_to_compare
    %find all rows in this condition, for this roitype, with non-nans in this datatype
    rows = find(all_long_vals(:,9,1)==c & all_long_vals(:,10,1)==4 & ~isnan(all_long_vals(:,1,1)));
    tmp = permute(all_long_vectors(:,1,:),[1 3 2]); %[measurement, vector]
    nr = length(rows);
    cur_mouse_vals = all_long_vals(rows,8,1); %store for long format
    cur_post_vals = squeeze(mean(tmp(rows,neuron_response_inds),2,'omitnan')); %delta mean dF/F during stim
    cur_pre_vals = squeeze(mean(tmp(rows-1,neuron_response_inds),2,'omitnan'));
    cur_delta_vals = cur_post_vals - cur_pre_vals;  %delta min dV/V during/after stim

    %store deltas
    if nr>size(delta_vals,1)
        pre_vals(end+1:nr,:,:) = nan;
        post_vals(end+1:nr,:,:) = nan;
        delta_vals(end+1:nr,:,:) = nan;
    end
    pre_vals(1:nr,c,1) = cur_pre_vals(:,:,1);
    post_vals(1:nr,c,1) = cur_post_vals(:,:,1);
    delta_vals(1:nr,c,1) = cur_delta_vals(:,:,1);

    %store measurement and metadata in long format for stats
    Data = [Data; cur_pre_vals];
    Mouse = [Mouse; cur_mouse_vals];
    Treatment = [Treatment; c*ones(nr,1)];
    Timepoint = [Timepoint; ones(nr,1)];
    Cell = [Cell; ones(nr,1)];

    Data = [Data; cur_post_vals];
    Mouse = [Mouse; cur_mouse_vals];
    Treatment = [Treatment; c*ones(nr,1)];
    Timepoint = [Timepoint; 2*ones(nr,1)];
    Cell = [Cell; ones(nr,1)];
end

%create figure for pre/post/delta boxplots
boxplot_vals = delta_vals(:,conds_to_compare([1 2]));
boxplot_colors = colors(conds_to_compare([1 2]),:);

if do_plots
    figure('Position',[100 100 box_plot_size])
    plot([0 6],[0 0],':k')
    hold on
    x = kron([1 2],ones(size(boxplot_vals,1),1));
    for j = 1:2
        swarmchart(x(:,j),boxplot_vals(:,x(1,j)),3,0.8*boxplot_colors(j,:),'filled','o','XJitterWidth',0.25,'MarkerFaceAlpha',0.7,'MarkerEdgeColor','none')
        b = boxchart(x(:,j), boxplot_vals(:,x(1,j)), 'BoxEdgeColor', boxplot_colors(j,:),'BoxFaceColor', 'none', 'LineWidth',0.5);
        b.MarkerStyle  = 'o';
        b.MarkerColor = SN_colors(5,:);
        b.MarkerSize = 0.5; 
    end

    plot([1 2],[box_ylimits(rt,2,dt) box_ylimits(rt,2,dt)],'k')
    plot([1 1],box_ylimits(rt,2,dt) - [0.05 0]*diff(box_ylimits(rt,:,dt)),'k')
    plot([2 2],box_ylimits(rt,2,dt) - [0.05 0]*diff(box_ylimits(rt,:,dt)),'k')
    % 
    plot([4 5],[box_ylimits(rt,2,dt) box_ylimits(rt,2,dt)],'k')
    plot([4 4],box_ylimits(rt,2,dt) - [0.05 0]*diff(box_ylimits(rt,:,dt)),'k')
    plot([5 5],box_ylimits(rt,2,dt) - [0.05 0]*diff(box_ylimits(rt,:,dt)),'k')

    text(1.2,box_ylimits(rt,2,dt) + 0.05*diff(box_ylimits(rt,:,dt)),'n.s. (p=0.00)','FontSize',4.0,'VerticalAlignment','bottom','HorizontalAlignment','left')
    text(4.2,box_ylimits(rt,2,dt) + 0.05*diff(box_ylimits(rt,:,dt)),'n.s. (p=0.00)','FontSize',4.0,'VerticalAlignment','bottom','HorizontalAlignment','left')
            
    ylim(box_ylimits(dt,:))
    xlim([0.25 2.75])
    set(gca,'FontSize',font_size)
    set(gca,'FontName','Arial')
    % ylabel(ylabels_boxplots{dt},'FontSize',font_size)
    ylabel('')
    xticks([1 2])
    xticklabels(xlabels_boxplots)
    ax = gca;
    % ax.Clipping = "off";
    box off

    print(gcf,'-vector','-dsvg',fullfile(base_folder,['2p boxplot ' dtypenames{dt} '.svg']))
end

if do_stats
    %run LME for this dataset
    cur_comparison = cur_comparison+1;
    % if dt==1 || dt==2 || dt==4 %neurons and vessel widths are 
    %     lmedata.measurement = log(Data(:,:,i) - min(Data(:,:,i)) + 0.1);
    % else
        lmedata.measurement = Data;
    % end
    lmedata.mouse = categorical(Mouse);
    lmedata.treatment = categorical(Treatment);
    lmedata.timepoint = categorical(Timepoint);
    lmedata.cell = categorical(Cell);
    lmetable = struct2table(lmedata);
    % treatment = fixed effect, mouse = random intercept
    % lme_model = fitlme(lmetable, 'measurement ~ treatment + (1|mouse)');
    lme_model = fitlme(lmetable, 'measurement ~ treatment*timepoint + (1|mouse) + (1|cell:mouse)');
    % ANOVA p-value for treatment effect
    anova_results = anova(lme_model); 
    lme_treatment_ps(dt,1,1) = lme_model.Coefficients.pValue(4); %post-psilocybin

    % disp([num2str(dt) 'AIC:']); 
    % disp(lme_model.ModelCriterion.AIC);

    %inspect lme model resiudals
    lme_residuals = residuals(lme_model);
    figure;
    subplot(2,2,1);
    histfit(lme_residuals); % Histogram with normal curve overlay
    title('residual distribution');
    subplot(2,2,2);
    qqplot(lme_residuals);
    title('Q-Q plot');
    subplot(2,2,3);
    plot(fitted(lme_model), lme_residuals, 'o');
    refline(0,0);
    xlabel('fitted values');
    ylabel('residuals');
    subplot(2,2,4);
    normplot(lme_residuals);
    title('normplot');
end


%%
%all other box plots and stats

%all_long_vectors:
%[measurement, datatype, vector]
%datatypes:
%1. dF/F
%2. Hb
%3. HbO
%4. dS/S
%5. cmro
%6. hbt
%7. width
%8. mouseID
%9. condition
%10. roitype (1=artery, 2=vein, 3=parenchyma, 4=gcamp-parenchyma)
%11. roinum

%set time ranges for analysis of neuron/vessel response
neuron_response_inds = find(t>0.5 & t<=6);
vessel_response_inds = find(t>0.5 & t<=12);
vessel_rise_inds = find(t>0.5 & t<=2.5);
vessel_recovery_inds = find(t>2.5 & t<=12);

do_stats = true;
do_plots = true;

box_plot_size = [120 60]; %smaller version
box_ylimits = nan(3,2,6); %[rt, x/y, dt]
%psilo
% box_ylimits(:,:,2) = [-15 15; -15 15; -10 10];
% box_ylimits(:,:,3) = [-10 10; -10 10; -5 5];
% box_ylimits(:,:,4) = [-20 20; -20 20; -10 10];
% box_ylimits(:,:,5) = [-20 20; -15 15; -10 10];
% box_ylimits(:,:,6) = [-25 25; -5 5; -10 10];
% conds_to_compare = [2 4];  %control pre/post, psilo pre/post
% xlabels_boxplots = {'ctrl','psil'};
%DOI
box_ylimits(:,:,2) = [-15 15; -15 15; -10 10];
box_ylimits(:,:,3) = [-12 12; -12 12; -7 7];
box_ylimits(:,:,4) = [-30 30; -20 20; -10 10];
box_ylimits(:,:,5) = [-21 21; -15 15; -10 10];
box_ylimits(:,:,6) = [-10 10; -10 10; -5 5];
conds_to_compare = [2 6];  %control pre/post, DOI pre/post
xlabels_boxplots = {'ctrl','DOI'};

cur_comparison = 0;
dtypenames = {'calcium','flow', 'hb','hbo', 'hbt','width'};
ylabels_boxplots = {'arterioles', 'venules', 'parenchyma'};
curlabel = 0;

cur_dt = 0;
for dt = [2 3 4 5 6] %calcium, flow, Hb, HbO, HbT, W
    num_rt = 3;
    if dt==6
        num_rt=2;
    end
    for rt = 1:num_rt
        pre_vals = nan(1,8,2);
        post_vals = nan(1,8,2);
        delta_vals = nan(1,8,2);
        Data = [];
        Mouse = [];
        Treatment = [];
        Timepoint = [];
        Cell = [];
        for c = conds_to_compare
            switch dt
                case 2
                    rows = find(all_long_vals(:,9,1)==c & all_long_vals(:,10,1)==rt & ~isnan(all_long_vals(:,4,1)));
                    tmp = permute(all_long_vectors(:,4,:),[1 3 2]); %[measurement, vector]
    
                case 3
                    rows = find(all_long_vals(:,9,1)==c & all_long_vals(:,10,1)==rt & ~isnan(all_long_vals(:,2,1)));
                    tmp = permute(all_long_vectors(:,2,:),[1 3 2])*1000000; %[measurement, vector]
    
                case 4
                    rows = find(all_long_vals(:,9,1)==c & all_long_vals(:,10,1)==rt & ~isnan(all_long_vals(:,3,1)));
                    tmp = permute(all_long_vectors(:,3,:),[1 3 2])*1000000; %[measurement, vector]

                case 5
                    rows = find(all_long_vals(:,9,1)==c & all_long_vals(:,10,1)==rt & ~isnan(all_long_vals(:,6,1)));
                    tmp = permute(all_long_vectors(:,6,:),[1 3 2])*1000000; %[measurement, vector]

                case 6
                    rows = find(all_long_vals(:,9,1)==c & all_long_vals(:,10,1)==rt & ~isnan(all_long_vals(:,7,1)));
                    tmp = permute(all_long_vectors(:,7,:),[1 3 2]); %[measurement, vector]
            end
            nr = length(rows);
            cur_mouse_vals = all_long_vals(rows,8,1); %store for long format
            cur_post_vals = squeeze(mean(tmp(rows,vessel_rise_inds),2,'omitnan')); 
            cur_pre_vals = squeeze(mean(tmp(rows-1,vessel_rise_inds),2,'omitnan'));
            cur_post_vals(:,:,2) = squeeze(mean(tmp(rows,vessel_recovery_inds),2,'omitnan'));
            cur_pre_vals(:,:,2) = squeeze(mean(tmp(rows-1,vessel_recovery_inds),2,'omitnan'));
            cur_delta_vals = cur_post_vals - cur_pre_vals; 

            %store deltas
            if nr>size(delta_vals,1)
                pre_vals(end+1:nr,:,:) = nan;
                post_vals(end+1:nr,:,:) = nan;
                delta_vals(end+1:nr,:,:) = nan;
            end
            pre_vals(1:nr,c,:) = cur_pre_vals;
            post_vals(1:nr,c,:) = cur_post_vals;
            delta_vals(1:nr,c,:) = cur_delta_vals;
    
            %store measurement and metadata in long format for stats
            Data = [Data; cur_pre_vals];
            Mouse = [Mouse; cur_mouse_vals];
            Treatment = [Treatment; c*ones(nr,1)];
            Timepoint = [Timepoint; ones(nr,1)];
            Cell = [Cell; all_long_vals(rows,11,1)];

            Data = [Data; cur_post_vals];
            Mouse = [Mouse; cur_mouse_vals];
            Treatment = [Treatment; c*ones(nr,1)];
            Timepoint = [Timepoint; 2*ones(nr,1)];
            Cell = [Cell; all_long_vals(rows,11,1)];
        end
        
        if do_stats
            for r_region = 1:2 %rise, recovery
                % lmedata.measurement = log(Data(:,:,i) - min(Data(:,:,i)) + 0.1);
                lmedata.measurement = Data(:,:,r_region);
                lmedata.mouse = categorical(Mouse);
                lmedata.treatment = categorical(Treatment);
                lmedata.timepoint = categorical(Timepoint);
                lmedata.cell = categorical(Cell);
                lmetable = struct2table(lmedata);
                % treatment = fixed effect, mouse = random intercept
                lme_model = fitlme(lmetable, 'measurement ~ treatment*timepoint + (1|mouse) + (1|cell:mouse)');
                % ANOVA p-value for treatment effect
                anova_results = anova(lme_model); 
                lme_treatment_ps(dt,rt,r_region) = lme_model.Coefficients.pValue(4);
    
                % disp([num2str(dt) 'AIC:']); 
                % disp(lme_model.ModelCriterion.AIC);
    
                %inspect lme model resiudals
                lme_residuals = residuals(lme_model);
                figure;
                subplot(2,2,1);
                histfit(lme_residuals); % Histogram with normal curve overlay
                title('residual distribution');
                subplot(2,2,2);
                qqplot(lme_residuals);
                title('Q-Q plot');
                subplot(2,2,3);
                plot(fitted(lme_model), lme_residuals, 'o');
                refline(0,0);
                xlabel('fitted values');
                ylabel('residuals');
                subplot(2,2,4);
                normplot(lme_residuals);
                title('normplot');
            end
        end

        boxplot_vals = [delta_vals(:,conds_to_compare([1 2]),1) delta_vals(:,conds_to_compare([1 2]),2)];
        boxplot_colors = colors(conds_to_compare([1 2 1 2]),:);
        fprintf(['dt ' num2str(dt) ' rt ' num2str(rt) ', min: ' num2str(min(boxplot_vals(:))), ' max: ' num2str(max(boxplot_vals(:))) '\n'])
        if do_plots
            figure('Position',[100 100 box_plot_size])
            plot([0 6],[0 0],':k')
            hold on
            x = kron([1 2 4 5],ones(size(boxplot_vals,1),1));
            for j = 1:4
                marker_size = 2;
                swarmchart(x(:,j),boxplot_vals(:,j),marker_size,0.6*boxplot_colors(j,:),'filled','o','XJitterWidth',0.25,'MarkerFaceAlpha',0.7,'MarkerEdgeColor','none')
                b = boxchart(x(:,j), boxplot_vals(:,j), 'BoxEdgeColor', boxplot_colors(j,:),'BoxFaceColor', 'none', 'LineWidth',0.5);
                b.MarkerStyle  = 'o';
                b.MarkerColor = SN_colors(5,:);
                b.MarkerSize = 0.5; 
            end
            plot([1 2],[box_ylimits(rt,2,dt) box_ylimits(rt,2,dt)],'k')
            plot([1 1],box_ylimits(rt,2,dt) - [0.05 0]*diff(box_ylimits(rt,:,dt)),'k')
            plot([2 2],box_ylimits(rt,2,dt) - [0.05 0]*diff(box_ylimits(rt,:,dt)),'k')
            % 
            plot([4 5],[box_ylimits(rt,2,dt) box_ylimits(rt,2,dt)],'k')
            plot([4 4],box_ylimits(rt,2,dt) - [0.05 0]*diff(box_ylimits(rt,:,dt)),'k')
            plot([5 5],box_ylimits(rt,2,dt) - [0.05 0]*diff(box_ylimits(rt,:,dt)),'k')
    
            text(1.2,box_ylimits(rt,2,dt) + 0.05*diff(box_ylimits(rt,:,dt)),'n.s. (p=0.00)','FontSize',4.0,'VerticalAlignment','bottom','HorizontalAlignment','left')
            text(4.2,box_ylimits(rt,2,dt) + 0.05*diff(box_ylimits(rt,:,dt)),'n.s. (p=0.00)','FontSize',4.0,'VerticalAlignment','bottom','HorizontalAlignment','left')
            
            ylim(box_ylimits(rt,:,dt))
            xlim([0 6])
            set(gca,'FontSize',font_size)
            set(gca,'FontName','Arial')
            % ylabel(ylabels_boxplots{rt},'FontSize',font_size)
            ylabel(' ','FontSize',font_size)
            xticks([1 2 4 5])
            % yticks(-50:10:50)
            xticklabels([xlabels_boxplots xlabels_boxplots])
            ax = gca;
            ax.Clipping = "off";
            box off
    
            print(gcf,'-vector','-dsvg',fullfile(base_folder,['2p boxplot ' ylabels_boxplots{rt} ' ' dtypenames{dt} '.svg']))
        end
    end
end

if do_stats
    % lme_treatment_ps
    % [~, ~, ~, adj_ps] = fdr_bh(lme_treatment_ps) %correct for all comparisons

    %remove nans before correction
    inds = nan(size(lme_treatment_ps));
    inds(1:end) = 1:numel(inds);
    tmpps = lme_treatment_ps;
    tmpps(isnan(lme_treatment_ps)) = [];
    tmpinds = inds;
    tmpinds(isnan(lme_treatment_ps)) = [];
    [~, ~, ~, tmp_adj_ps] = fdr_bh(lme_treatment_ps); %correct for all comparisons
    adj_ps = nan(size(lme_treatment_ps));
    adj_ps(tmpinds) = tmpps
end



%% plot baseline vessel widths (post/pre)
%all_long_vectors:
%[measurement, datatype, vector]
%datatypes:
%1. dF/F
%2. Hb
%3. HbO
%4. dS/S
%5. cmro
%6. hbt
%7. width
%8. mouseID
%9. condition
%10. roitype (1=artery, 2=vein, 3=parenchyma, 4=gcamp-parenchyma)
%11. roinum

upp = 2.68; %microns/pixel
box_plot_size = [160 80];
baseline_ylimits = [-50 65];
lme_baseline_ps = nan(5,1);
dtypenames = {'calcium','flow', 'hb','hbo', 'hbt','width'};
ylabels_boxplots = {'\Delta response (% \DeltaF/F)', '\Delta recovery (% \DeltaS/S)', '\Delta recovery (\DeltaHb)', '\Delta recovery (\DeltaHbO)', '\Delta recovery (\DeltaHbT)', 'baseline % \DeltaW/W'};
curlabel = 0;
boxplot_vals = nan(1,6);

do_stats = false;
% conds_to_compare = [2 4];  %control pre/post, psilo pre/post
% conds_to_compare = [2 6];  %control pre/post, DOI pre/post
% conds_to_compare = [4 6];  %psilo pre/post, DOI pre/post

do_plots = true;
conds_to_compare = [2 4 6];  %control pre/post, psilo pre/post, DOI pre/post

dt = 6;
for rt = 1:2
    pre_vals = nan(1,8);
    post_vals = nan(1,8);
    delta_vals = nan(1,8);
    Data = [];
    Treatment = [];
    Mouse = [];
    Timepoint = [];
    Cell = [];
    
    for c = conds_to_compare
        rows = find(all_long_vals(:,9,1)==c & all_long_vals(:,10,1)==rt & ~isnan(all_long_vals(:,12,1)));
        tmp = permute(all_long_vectors(:,12,:),[1 3 2]); %[measurement, vector]
        nr = length(rows);
        cur_mouse_vals = all_long_vals(rows,8,1); %store for long format
        cur_post_vals = squeeze(mean(tmp(rows,baseline_frames),2,'omitnan'));
        cur_pre_vals = squeeze(mean(tmp(rows-1,baseline_frames),2,'omitnan'));
        cur_delta_vals = 100*((cur_post_vals-cur_pre_vals)./cur_pre_vals);

        %store deltas
        % cur_delta_vals
        if nr>size(delta_vals,1)
            pre_vals(end+1:nr,:,:) = nan;
            post_vals(end+1:nr,:,:) = nan;
            delta_vals(end+1:nr,:,:) = nan;
        end
        pre_vals(1:nr,c) = cur_pre_vals;
        post_vals(1:nr,c) = cur_post_vals;
        delta_vals(1:nr,c) = cur_delta_vals;

        %store measurement and metadata in long format for stats
        Data = [Data; cur_pre_vals];
        Mouse = [Mouse; cur_mouse_vals];
        Treatment = [Treatment; c*ones(nr,1)];
        Timepoint = [Timepoint; ones(nr,1)];
        Cell = [Cell; all_long_vals(rows,11,1)];

        Data = [Data; cur_post_vals];
        Mouse = [Mouse; cur_mouse_vals];
        Treatment = [Treatment; c*ones(nr,1)];
        Timepoint = [Timepoint; 2*ones(nr,1)];
        Cell = [Cell; all_long_vals(rows,11,1)];
    end

    if do_stats
        %run LME for this dataset
        cur_comparison = cur_comparison+1;
        % lmedata.measurement = log(Data - min(Data) + 0.01);
        lmedata.measurement = Data;
        lmedata.mouse = categorical(Mouse);
        lmedata.treatment = categorical(Treatment);
        lmedata.timepoint = categorical(Timepoint);
        lmedata.cell = categorical(Cell);
        lmetable = struct2table(lmedata);
        % treatment = fixed effect, mouse = random intercept
        % lme_model = fitlme(lmetable, 'measurement ~ treatment + (1|mouse)');
        lme_model = fitlme(lmetable, 'measurement ~ treatment*timepoint + (1|mouse) + (1|cell:mouse)');
                
        % ANOVA p-value for treatment effect
        % anova_results = anova(lme_model); 
        lme_baseline_ps(rt,:) = lme_model.Coefficients.pValue(4);

        % disp([num2str(dt) 'AIC:']); 
        % disp(lme_model.ModelCriterion.AIC);
        

        %inspect lme model resiudals
        lme_residuals = residuals(lme_model);
        % disp(skewness(residuals(lme_model)));
        figure;
        subplot(2,2,1);
        histfit(lme_residuals); % Histogram with normal curve overlay
        title('residual distribution');
        subplot(2,2,2);
        qqplot(lme_residuals);
        title('Q-Q plot');
        subplot(2,2,3);
        plot(fitted(lme_model), lme_residuals, 'o');
        refline(0,0);
        xlabel('fitted values');
        ylabel('residuals');
        subplot(2,2,4);
        normplot(lme_residuals);
        title('normplot');
    end

    if do_plots
        nr = size(delta_vals,1);
        if nr>size(boxplot_vals,1)
            boxplot_vals(end+1:nr,:,1) = nan;
        end
        if rt==1
            boxplot_vals(1:nr,1:3) = delta_vals(:,conds_to_compare([1 2 3]),1);
        end
        if rt==2
            boxplot_vals(1:nr,4:6) = delta_vals(:,conds_to_compare([1 2 3]),1);
        end
    end
end

if do_plots
    boxplot_colors = colors(conds_to_compare([1 2 3 1 2 3]),:);

    figure('Position',[100 100 box_plot_size])
    plot([0 9],[0 0],':k')
    hold on
    x = kron([1 2 3 5 6 7],ones(size(boxplot_vals,1),1));

    for j = 1:6
        marker_size = 3;
        swarmchart(x(:,j),boxplot_vals(:,j),marker_size,0.6*boxplot_colors(j,:),'filled','o','XJitterWidth',0.25,'MarkerFaceAlpha',0.7,'MarkerEdgeColor','none')
        b = boxchart(x(:,j), boxplot_vals(:,j), 'BoxEdgeColor', boxplot_colors(j,:),'BoxFaceColor', 'none', 'LineWidth',0.5);
        b.MarkerStyle  = 'o';
        b.MarkerColor = SN_colors(5,:);
        b.MarkerSize = marker_size/4; 
    end

    ylim(baseline_ylimits)
    xlim([0 8])
    set(gca,'FontSize',font_size)
    set(gca,'FontName','Arial')
    ylabel(ylabels_boxplots{dt},'FontSize',font_size)
    xticks([1 2 3 5 6 7])
    % yticks(-50:10:50)
    xticklabels({'ctrl','psil','DOI','ctrl','psil','DOI'})
    % xticklabels({'p','p','d','','p','p','d'})
    ax = gca;
    ax.Clipping = "off";
    box off

    print(gcf,'-vector','-dsvg',fullfile(base_folder,['2p baseline boxplot ' dtypenames{dt} '.svg']))
end

%%
lme_baseline_ps(1) = 0.2293; %control vs psilo, Widefield arterioles
lme_baseline_ps(2) = 0.0003; %control vs psilo, Widefield venules
lme_baseline_ps(3) = 0.0005156; %control vs DOI, Widefield arterioles
lme_baseline_ps(4) = 0.0000; %control vs DOI, Widefield venules
lme_baseline_ps(5) = 0.0246; %psilo vs DOI, Widefield arterioles
lme_baseline_ps(6) = 0.0041; %psilo vs DOI, Widefield venules
lme_baseline_ps(7) = 0.0031562; %control vs psilo, from 2P capillaries
lme_baseline_ps(8) = 0.36473; %control vs psilo+MDL, from 2P capillaries
lme_baseline_ps(9) = 0.082274; %psilo vs psilo+MDL, from 2P capillaries

% lme_baseline_ps
[~, ~, ~, adj_b_ps] = fdr_bh(lme_baseline_ps)



%% create image of baseline width measurement

I = read_file(fullfile(data_folder,'raw_ref.tif'));
%speckle, gcamp, 530, 560
image_size = [68 68];
i = 3;
tmp = I(:,:,i);
tmp = tmp-prctile(tmp(:),2);
tmp = tmp/prctile(tmp(:),98);
if i>2
    tmp = tmp*2.5;
end
figure('Position',[100 100 image_size])
ax = axes('Parent', gcf, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
imshow(tmp)
hold on
plot(boundary_x,boundary_y,'--w','LineWidth',1)

% create patches of artery ROIs
load(fullfile(data_folder, 'vesselWidths_vessels_ref.mat'));
for v = 1:vesselBoxes.num
    x = vesselBoxes.x(v,:);
    y = vesselBoxes.y(v,:);
    switch vesselBoxes.type{v}
        case 'Artery'
            patch(x,y,'k','FaceColor',SN_colors(5,:),'EdgeColor',SN_colors(5,:),'FaceAlpha',0,'LineWidth',1)
        case 'Vein'
            patch(x,y,'k','FaceColor',SN_colors(10,:),'EdgeColor',SN_colors(10,:),'FaceAlpha',0,'LineWidth',1)
    end
end

print(gcf,'-vector','-dsvg',fullfile(base_folder,['widefield image 1 ' image1names{i} '.svg']))

