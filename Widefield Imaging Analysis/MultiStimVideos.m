
%ask user for mat filename, get corresponding tif filename
if exist('path','var') && ischar(path)
    [file,path] = uigetfile([path '*.mat']);
else
    [file,path] = uigetfile('*.mat');
end
matfile = [path file];
imfile = [matfile(1:find(matfile=='.',1,'last')) 'tif'];
savedir = fullfile(path,['MSV4_' file(1:end-4)]);
if ~exist(savedir,'dir')
    mkdir(savedir)
end

TifLink = Tiff(imfile,'r');
load(matfile,'trigTime','FrameData','wlopts','SessionInfo','Xq','trigDuration');
%Xq: 10 Hz time samples per stimulus (-3 to 16.6 s)
% wlopts: {'Speckle', 'Fluorescence', 'Reflectance', 'Reflectance'}
% WavelengthNumbers: [1=785, 2=470; 3=530; 4=560]
%trigTime: timepoints of triggers (in s)
%trigDuration: length of trigger sections (in s)
output_names = {'speckle','gcamp','hb','hbo','cmro','raw530','cbf'};
sz = [SessionInfo.ROI([4 3])/SessionInfo.Binning length(Xq)]; %[900 1000 197]
refIdx = find(contains(wlopts,'Speckle'));
speckleExposure = SessionInfo.Exposure(refIdx);
Ntrig = length(trigTime); %number of triggers (15)
WNum = length(wlopts); %number of excitation wavelengths (4)

filter_size = 3;
filter_type = 'gaussian'; % 'gaussian' (size=sigma) or 'square' (flat; size=width)

%create Tiff tagstruct
tagstruct.ImageLength = sz(1);
tagstruct.ImageWidth = sz(2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Software = 'MATLAB';


%(part 1) load and process Speckle, Hb, HbO, and GCaMP data 
%prepare to write data to .tiff stacks
TifWrite_sp = Tiff(fullfile(savedir,'speckle.tif'),'w8');
TifWrite_gc = Tiff(fullfile(savedir,'gcamp.tif'),'w8');
TifWrite_hb = Tiff(fullfile(savedir,'hb.tif'),'w8');
TifWrite_hbo = Tiff(fullfile(savedir,'hbo.tif'),'w8');
TifWrite_cmro = Tiff(fullfile(savedir,'cmro.tif'),'w8');
TifWrite_raw530 = Tiff(fullfile(savedir,'raw530.tif'),'w8');
TifWrite_cbf = Tiff(fullfile(savedir,'cbf.tif'),'w8');
tifstart = true;

f = waitbar(0,'Processing full speckle, Hb, HbO, GCaMP, CMRO, CBF videos (step 1/5)');
for tn=1:Ntrig
    waitbar(0.5*(tn/Ntrig),f,'Processing full speckle, Hb, HbO, GCaMP, CMRO videos (step 1/5)');

    % get baselines for all wavelengths
    Baselines = zeros([sz(1:2) 5]);
    for w = 1:4 %[785, 470, 530, 560]
        %get baseline image by averaging all baseline frames
        baseidx = find(FrameData.Time>=trigTime(tn)-3.2 & FrameData.Time<=trigTime(tn)+0.2 & FrameData.WavelengthNumber==w);
        baseL = length(baseidx);
        for i=1:baseL
            TifLink.setDirectory(baseidx(i));
            if w==1
                tmp_frame = Pipeline.FastSpeckle(TifLink.read); %read speckle frame
                Baselines(:,:,1) = Baselines(:,:,1) + tmp_frame; %add baseline frame for speckle
                Baselines(:,:,5) = Baselines(:,:,5) + 1./tmp_frame.^2.*(1e6./speckleExposure); %add baseline frame for ICT
            else
                Baselines(:,:,w) = Baselines(:,:,w) + double(TifLink.read); %add frame
            end
        end
        Baselines(:,:,w) = Baselines(:,:,w)/baseL; %divide by number of add frames to get mean
        switch filter_type
            case 'gaussian'
                Baselines(:,:,w) = imgaussfilt(Baselines(:,:,w),filter_size,'padding','replicate');
            case 'square'
                Baselines(:,:,w) = conv2(Baselines(:,:,w),(1/(filter_size^2))*ones(filter_size),'same');
            otherwise
                error(['unknown filter_type "' filter_type '"'])
        end
    end

    %get data relative to baselines for all wavelengths for the entire triggered region
    aroundFrames = nan([sz(1:2) 4 4]); %[H, W, frame, wavelength] actual frames before and after current time
    preFrames = nan([sz(1:2) 4]); % temporally-smoothed frame BEFORE current time
    postFrames = nan([sz(1:2) 4]); % temporally-smoothed frame AFTER current time
    curFrames = nan([sz(1:2) 4]); % interpolated frame AT current time
    aroundFramesidx = nan(4,4); %[wavelength, frame]
    previous_aroundFramesidx = zeros(4,4); %[wavelength, frame]
    for i = 1:sz(3)
        curTime = Xq(i);

        %get frame for current timepoint (interpolate between nearest, temporally-smoothed frames)
        for w = 1:4
            %find frames for this wavelength that surround the current timepoint (2 before and 2 after)
            aroundFramesidx(w,:) = [find(FrameData.Time<=trigTime(tn)+curTime & FrameData.WavelengthNumber==w, 2, "last"); ...
                find(FrameData.Time>trigTime(tn)+curTime & FrameData.WavelengthNumber==w, 2, "first")];
            
            %check if any of these frames have already been loaded
            if all(aroundFramesidx(w,:)==previous_aroundFramesidx(w,:))
                %aroundFramees is unchanged, therefore pre/postFrames are unchanged
            elseif all(aroundFramesidx(w,1:3)==previous_aroundFramesidx(w,2:4))
                %shifted to one new frame; old postFrames is new preFrame, but new preFrame needs to be calculated
                preFrames(:,:,w) = postFrames(:,:,w);
                aroundFrames(:,:,1:3,w) = aroundFrames(:,:,2:4,w);
                TifLink.setDirectory(aroundFramesidx(w,4));
                if w==1
                    aroundFrames(:,:,4,w) = Pipeline.FastSpeckle(TifLink.read);
                else
                    aroundFrames(:,:,4,w) = double(TifLink.read);
                end
                postFrames(:,:,w) = mean(aroundFrames(:,:,2:4,w),3,'omitnan');
            else
                %multiple new frames need to be loaded; both pre and post frames need to be recalculated
                for a = 1:4
                    if any(aroundFramesidx(w,a)==previous_aroundFramesidx(w,:))
                        %keep already loaded frame
                        old_a = find(previous_aroundFramesidx(w,:)==aroundFramesidx(w,a));
                        aroundFrames(:,:,a,w) = aroundFrames(:,:,old_a,w);
                    else
                        %load new frame
                        TifLink.setDirectory(aroundFramesidx(w,a));
                        if w==1
                            aroundFrames(:,:,a,w) = Pipeline.FastSpeckle(TifLink.read);
                        else
                            aroundFrames(:,:,a,w) = double(TifLink.read);
                        end
                    end
                end
                preFrames(:,:,w) = mean(aroundFrames(:,:,1:3,w),3,'omitnan');
                postFrames(:,:,w) = mean(aroundFrames(:,:,2:4,w),3,'omitnan');
            end
            
            %linearly interpolate between pre and post frames to get current frame
            preWeight = (FrameData.Time(aroundFramesidx(w,3)) - (trigTime(tn)+curTime))/(FrameData.Time(aroundFramesidx(w,3)) - FrameData.Time(aroundFramesidx(w,2)));
            curFrames(:,:,w) = preWeight*postFrames(:,:,w) + (1-preWeight)*postFrames(:,:,w);
            if w==3
                raw530 = curFrames(:,:,w);
            end
            switch filter_type
                case 'gaussian'
                    curFrames(:,:,w) = imgaussfilt(curFrames(:,:,w),filter_size,'padding','replicate');
                case 'square'
                    curFrames(:,:,w) = conv2(curFrames(:,:,w),(1/(filter_size^2))*ones(filter_size),'same');
                otherwise
                    error(['unknown filter_type "' filter_type '"'])
            end
        end

        %get speckle contrast (wavelength 1)
        ImageSp = (Baselines(:,:,1)./curFrames(:,:,1)).^2;

        %get log(I0/I) for 530, 560 (wavelengths 3-4)
        Reflectance530 = curFrames(:,:,3)./Baselines(:,:,3);
        Image530 = log(Baselines(:,:,3)./curFrames(:,:,3));
        Image560 = log(Baselines(:,:,4)./curFrames(:,:,4));

        %extract Hb, HbO from 530/560 log(I0/I) data
        coeff = [41778 39735;...   % 530 [HbO Hb]
                  32617 53817];     % 560 [HbO Hb]
        dpf = [0.0371713; 0.0392168]; % expected pathlength factors in cm for proper conversion to uM
        cx = coeff.*dpf;
        invM = cx^-1;
        ImageHbO = invM(1,1)*Image530 + invM(1,2)*Image560;
        ImageHb = invM(2,1)*Image530 + invM(2,2)*Image560;

        %get normalized gcamp fluorescence (wavelength 2)
        ImageGc = curFrames(:,:,2)./Baselines(:,:,2);

        %correct normalized gcamp fluorescence with HbO
        Coeffs = [33767 16450;...   % 470
                  41778 39735];     % 530
        Xex = 0.056; Xem = 0.057;
        ua_ex = Coeffs(1,1)*ImageHbO+Coeffs(1,2)*ImageHb;
        ua_em = Coeffs(2,1)*ImageHbO+Coeffs(2,2)*ImageHb;
        ImageGc = ImageGc.*exp(ua_ex*Xex+ua_em*Xem);

        % alternatively use reflectance for simpler gcamp correction
        % ImageGc = ImageGc./Reflectance530;

        %get ICT, dCBF, CMRO2
        imageICT = 1./curFrames(:,:,1).^2.*(1e6./speckleExposure);
        imageCBF = imageICT./Baselines(:,:,5);
        HbT0 = 100e-6;
        HbR0 = 60e-6;
        imageCMRO = 1./imageCBF.*(1+(ImageHb+ImageHbO)/HbT0).*(1+ImageHb./HbR0)-1;
        
        %save current frame indices for next loop
        previous_aroundFramesidx = aroundFramesidx;
        
        %append all images (Speckle, Hb, HbO, GCaMP) to tif stacks
        if tifstart
            tifstart = false;            
        else
            TifWrite_sp.writeDirectory();
            TifWrite_hb.writeDirectory();
            TifWrite_hbo.writeDirectory();
            TifWrite_gc.writeDirectory();
            TifWrite_cmro.writeDirectory();
            TifWrite_raw530.writeDirectory();
            TifWrite_cbf.writeDirectory();
        end
        TifWrite_sp.setTag(tagstruct);
        TifWrite_hb.setTag(tagstruct);
        TifWrite_hbo.setTag(tagstruct);
        TifWrite_gc.setTag(tagstruct);
        TifWrite_cmro.setTag(tagstruct);
        TifWrite_raw530.setTag(tagstruct);
        TifWrite_cbf.setTag(tagstruct);
        TifWrite_sp.write(single(ImageSp));
        TifWrite_hb.write(single(ImageHb));
        TifWrite_hbo.write(single(ImageHbO));
        TifWrite_gc.write(single(ImageGc));
        TifWrite_cmro.write(single(imageCMRO));
        TifWrite_raw530.write(single(raw530));
        TifWrite_cbf.write(single(imageCBF));
    end
end
TifLink.close;
TifWrite_sp.close();
TifWrite_hb.close();
TifWrite_hbo.close();
TifWrite_gc.close();
TifWrite_cmro.close();
TifWrite_raw530.close();
TifWrite_cbf.close();
pause(1)
clear TifLink TifWrite_sp TifWrite_hb TifWrite_hbo TifWrite_gc TifWrite_cmro TifWrite_raw530 TifWrite_cbf
pause(2)



%(part 2) calculate average-trial videos of all full videos
for w = 1:length(output_names)
    waitbar(0.5 + 0.2*(w/5),f,'Calculating trial-averaged videos (step 2/5)');
    TifLink = Tiff(fullfile(savedir,[output_names{w} '.tif']),'r');
    TifWrite = Tiff(fullfile(savedir,[output_names{w} '_mean.tif']),'w8');
    tifstart = true;
    for i = 1:sz(3)
        CurImage = zeros(sz(1:2));
        for t = 1:Ntrig
            %load frame i for trig t
            TifLink.setDirectory(i + sz(3)*(t-1));
            CurImage = CurImage + double(TifLink.read);
        end
    
        %save mean frame i for all trigs
        CurImage = CurImage/Ntrig;
        if tifstart
            tifstart = false;            
        else
            TifWrite.writeDirectory();
        end
        TifWrite.setTag(tagstruct);
        TifWrite.write(single(CurImage));
    end
    TifLink.close;
    TifWrite.close();
end



% (part 3) delete full videos
for w = 1:length(output_names)
    waitbar(0.7+0.1*(w/5),f,'Deleting full videos (step 3/5)');
    filename = fullfile(savedir,[output_names{w} '.tif']);
    delete(filename)
end



% (part 4) create mean images of all 4 raw channels
tifstart = true;
TifWrite = Tiff(fullfile(savedir,'raw_ref.tif'),'w8');
for w = 1:4 %[785, 470, 530, 560]
    waitbar(0.8+0.1*(w/4),f,'Creating mean images of raw imaging channels (step 4/5)');
    %get random subset of 100 frames 
    w_inds = find(FrameData.WavelengthNumber==w);
    wL = min([200 length(w_inds)]);
    w_inds = w_inds(randperm(length(w_inds),wL)); %get randomized subset (or whole set if <=200)
    CurImage = zeros(sz(1:2));
    TifLink = Tiff(imfile,'r');
    for i=1:wL
        TifLink.setDirectory(w_inds(i));
        CurImage = CurImage + double(TifLink.read); %add frame
    end
    CurImage = CurImage/wL; %divide by number of add frames to get mean
    if tifstart
        tifstart = false;            
    else
        TifWrite.writeDirectory();
    end
    TifWrite.setTag(tagstruct);
    TifWrite.write(single(CurImage));
    TifLink.close;

    if w==3 %save 530 nm separately (nice view of vessels)
        TifWrite2 = Tiff(fullfile(savedir, 'vessels_ref.tif'),'w8');
        TifWrite2.setTag(tagstruct);
        TifWrite2.write(single(CurImage));
        TifWrite2.close();
    end
end
TifWrite.close();



% (part 5) create oxygenation map
%set timepoints for max Hb response (indicates veins)
max_response_inds = find(Xq>=1.5 & Xq<=5.5);
responseN = length(max_response_inds);

%get timepoints for early HbO response (indicates arteries)
tmp = find(Xq>=0.5 & Xq<=1.5);
responseN(2) = length(tmp);
max_response_inds(2,1:responseN(2)) = tmp;

for w = 3:4 %Hb, HbO
    waitbar(0.9+0.1*((w-2)/2),f,'Creating 2-channel oxygenation map (step 5/5)');
    TifLink = Tiff(fullfile(savedir,[output_names{w} '_mean.tif']),'r');
    CurImage = zeros(sz(1:2));
    for i = 1:responseN(w-2)
        ind = max_response_inds(w-2,i);
        TifLink.setDirectory(ind);
        CurImage = CurImage + double(TifLink.read);
    end
    CurImage = CurImage/responseN(w-2);

    if w==3 %flip, scale, and add Hb image to red channel (veins)
        CurImage = CurImage-prctile(CurImage,1,'all');
        CurImage = CurImage/prctile(CurImage,99,'all');
        CurImage = 1-CurImage;
        oxyImage = zeros([size(CurImage) 3]);
        oxyImage(:,:,1) = CurImage;
    else %scale and add HbO image to blue channel (arteries)
        CurImage = CurImage-prctile(CurImage,1,'all');
        CurImage = CurImage/prctile(CurImage,99,'all');
        oxyImage(:,:,3) = CurImage;
    end
    TifLink.close;
end
TifWrite = Tiff(fullfile(savedir, 'oxymap_ref.tif'),'w8');
%create Tiff tagstruct
tagstruct.ImageLength = sz(1);
tagstruct.ImageWidth = sz(2);
tagstruct.Photometric = Tiff.Photometric.RGB;
tagstruct.BitsPerSample = 32;
tagstruct.SamplesPerPixel = 3;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Software = 'MATLAB';
TifWrite.setTag(tagstruct);
TifWrite.write(single(oxyImage));
TifWrite.close();



%close waitbar
delete(f)