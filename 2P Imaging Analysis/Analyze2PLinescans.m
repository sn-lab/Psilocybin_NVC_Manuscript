
%% set base directory
base_folder = 'D:\2P Data\New file organization\';


%% (STEP 1) get data directories (new file organization)

% get velocity directories
cond_names = {'before-control', 'after-control', 'before-psilocybin', 'after-psilocybin', 'before-mdl', 'after-mdl', 'before-mdl-psilocybin', 'after-mdl-psilocybin'};
treatment_names = {'control', 'control', 'psilocybin', 'psilocybin', 'mdl', 'mdl', 'mdl-psilocybin', 'mdl-psilocybin'};
timepoint_names = {'pre','post'};
vel_roi_nums = [1 5 9]; %linescan ROI #s associated with down-vessel scans
width_roi_nums = [4 8 12]; %linescan ROI #s associated with cross-vessel scans
neuron_roi_nums = [2 3 6 7 10 11]; %linescan ROI #s associated with cross-neuron cell body scans

%set time vector variables
P_line = 0.006; %linescan period (s)
P = P_line*50; %period between sets of 50 linescans (size of blocks for radon transformation)
time_vec = -10*P:P:66*P; %~ -3 s to +19.8 s
baseline_inds = find(time_vec<0);
stim_inds = find(time_vec>0.5 & time_vec<6);
post_stim_inds = find(time_vec>6.5 & time_vec<12);
stim_and_post_stim_inds = find(time_vec>0.5 & time_vec<12);
num_pts = length(time_vec);

%search the base directory for all mouse data folders
l = dir(base_folder);
for i = length(l):-1:1
    %exclude non-mouse folders (mouse folders start with a number)
    if isnan(str2double(l(i).name(1)))
        l(i) = [];
    end
end
num_mice = length(l);
num_conds = length(cond_names);
mouse_names = cell(num_mice,1);
for i = 1:num_mice
    mouse_names{i} = l(i).name;
end

%search through data folders for mouse/condition pairs which contain flow data
num_vessels = length(vel_roi_nums);
velnames = cell(num_conds,num_mice,num_vessels);
foundvels = false(num_conds,num_mice,num_vessels);
for m = 1:num_mice
    for c = 1:num_conds
        cur_dir = fullfile(base_folder,l(m).name,cond_names{c},'stim');
        for v = 1:num_vessels
            cur_name = ls(fullfile(cur_dir,['*ROI' num2str(vel_roi_nums(v)) '*rawVel 50100.mat']));
            if ~isempty(cur_name)
                velnames{c,m,v} = fullfile(cur_dir,cur_name);
                foundvels(c,m,v) = true;
            end
        end
    end
end

%remove ROIs which weren't measured in both the before and after conditions
tmp_before = foundvels([1 3 5 7],:,:); %just the "before" conditions
tmp_after = foundvels([2 4 6 8],:,:); %just the "after" conditions
not_before_idx = tmp_before==0; %missing ROIs in the "before" conditions
not_after_idx = tmp_after==0; %missing ROIs in the "after" conditions
tmp_after(not_before_idx) = 0; %remove missing "befores" in the "afters"
tmp_before(not_after_idx) = 0; %remove missing "afters" in the "befores"
foundvels([1 3 5 7],:,:) = tmp_before;
foundvels([2 4 6 8],:,:) = tmp_after;

%search through data folders for mouse/condition pairs which contain neuron data
num_neurons = length(neuron_roi_nums);
neuronnames = cell(num_conds,num_mice,num_neurons);
foundneurons = false(num_conds,num_mice,num_neurons);
for c = 1:num_conds
    for m = 1:num_mice
        tmp = velnames{c,m,1};
        if ~isempty(tmp)
            [basedir,fname] = fileparts(tmp);
            for n = 1:num_neurons
                tmp = fullfile(basedir,strcat("ALNA1_ROI_",num2str(neuron_roi_nums(n)),"*.tif"));
                tmp = ls(tmp);
                if ~isempty(tmp)
                    neuronnames{c,m,n} = fullfile(basedir,tmp); 
                    foundneurons(c,m,n) = true;
                end
            end
        end
    end
end

%remove ROIs which weren't measured in both the before and after conditions
tmp_before = foundneurons([1 3 5 7],:,:); %just the "before" conditions
tmp_after = foundneurons([2 4 6 8],:,:); %just the "after" conditions
not_before_idx = tmp_before==0; %missing ROIs in the "before" conditions
not_after_idx = tmp_after==0; %missing ROIs in the "after" conditions
tmp_after(not_before_idx) = 0; %remove missing "befores" in the "afters"
tmp_before(not_after_idx) = 0; %remove missing "afters" in the "befores"
foundneurons([1 3 5 7],:,:) = tmp_before;
foundneurons([2 4 6 8],:,:) = tmp_after;


%search through data folders for mouse/condition pairs which contain vessel width data
diamnames = cell(num_conds,num_mice,num_vessels);
founddiams = false(num_conds,num_mice,num_vessels);
for c = 1:num_conds
    for m = 1:num_mice
        tmp = velnames{c,m,1};
        if ~isempty(tmp)
            [basedir,fname] = fileparts(tmp);
            for v = 1:num_vessels
                tmp = fullfile(basedir,strcat("AL2_Results_ROI_",num2str(width_roi_nums(v)),"*.mat"));
                tmp = ls(tmp);
                if ~isempty(tmp)
                    diamnames{c,m,v} = fullfile(basedir,tmp); 
                    founddiams(c,m,v) = true;
                end
            end
        end
    end
end

%remove ROIs which weren't measured in both the before and after conditions
tmp_before = founddiams([1 3 5 7],:,:); %just the "before" conditions
tmp_after = founddiams([2 4 6 8],:,:); %just the "after" conditions
not_before_idx = tmp_before==0; %missing ROIs in the "before" conditions
not_after_idx = tmp_after==0; %missing ROIs in the "after" conditions
tmp_after(not_before_idx) = 0; %remove missing "befores" in the "afters"
tmp_before(not_after_idx) = 0; %remove missing "afters" in the "befores"
founddiams([1 3 5 7],:,:) = tmp_before;
founddiams([2 4 6 8],:,:) = tmp_after;

fprintf(['Files found for ' num2str(num_mice) ' mice\n'])

%list of mice:
% 1-1L
% 2-2114
% 3-2117
% 4-2119
% 5-2236
% 6-2621
% 7-2R
% 8-3108
% 9-3109
% 10-3113
% 11-3114
% 12-3115
% 13-3815
% 14-3817
% 15-3904
% 16-3905
% 17-3N
% 18-4B
% 19-543959
% 20-546842
% 21-547880
% 22-6184N_M
% 23-6185R_M
% 24-6186
% 25-6187_F
% 26-6189
% 27-6680_None
% 28-6681_Top
% 29-6682_Bottom
% 30-7675
% 31-7676
% 32-885


%% (OPTIONAL) check for missing datasets
%some neurons/vessel measurements are missing for various reasons (e.g.
%poor data quality; could not find enough close neuron cell bodies)

foundsomething = sum(founddiams,3,'omitnan')>0 | sum(foundvels,3,'omitnan')>0 | sum(foundneurons,3,'omitnan')>0;
foundall = (sum(founddiams,3,'omitnan')>0) & (sum(foundvels,3,'omitnan')>0) & (sum(foundneurons,3,'omitnan')>0);
missingsomething = foundsomething & ~foundall;

for m = 1:num_mice
    for c = 1:num_conds
        if missingsomething(c,m)
            fprintf(['Mouse ' l(m).name ' cond ' num2str(c) ', missing '])
            if sum(foundvels(c,m,:),3,'omitnan')==0
                fprintf(' - vels ')
            end
            if sum(founddiams(c,m,:),3,'omitnan')==0
                fprintf(' - diams ')
            end
            if sum(foundneurons(c,m,:),3,'omitnan')==0
                fprintf(' - neurons ')
            end
            fprintf('\n')
        end
    end
end


%% (OPTIONAL) check for missing stimulus trigger files
missing_files = nan(1,2);
num_missing_files = 0;
for m = 1:num_mice
    for c = 1:num_conds
        if any([squeeze(foundvels(c,m,:)); squeeze(founddiams(c,m,:)); squeeze(foundneurons(c,m,:))])
            cur_dir = fullfile(base_folder,l(m).name,cond_names{c},'stim');
            curtrigfilename = ls(fullfile(cur_dir,'Ch1_triggers_extracted*'));
            if isempty(curtrigfilename)
                fprintf(['missing trigger file: mouse ' num2str(m) ', cond ' num2str(c) '\n'])
                num_missing_files = num_missing_files+1;
                missing_files(num_missing_files,1) = m;
                missing_files(num_missing_files,2) = c;
            end
        end
    end
end
if num_missing_files==0
    fprintf('All trigger files found\n')
else
    fprintf([num2str(num_missing_files) ' missing trigger files\n'])
end


%% (ALREADY DONE) create stimulus trigger files
numChannels = 4;
samplesPerFrame = 15000;
% samplesPerFrame = 14875; %for mouse 6681 (mistaken duration for 1 ROI)
for f = 1:num_missing_files
    m = missing_files(f,1);
    c = missing_files(f,2);
    fprintf('file %d of %d (m%d,c%d)\n',f,num_missing_files,m,c)

    %load pmt data
    cur_dir = fullfile(base_folder,l(m).name,cond_names{c},'stim');
    pmt_name = fullfile(cur_dir,ls(fullfile(cur_dir,'*.pmt.dat')));
    assert(~isempty(pmt_name),'cannot find pmt.dat file')
    fileID = fopen(pmt_name);

    %read data in chunks of numChannels*samplesPerFrame*1000 samples at a time
    count = numChannels*samplesPerFrame*linesToRead;
    linesToRead = 1000;
    channelAvg = [];
    while count==(numChannels*samplesPerFrame*linesToRead)
        [data, count] = fread(fileID,numChannels*samplesPerFrame*linesToRead,'int16');
    
        %extract channel 1
        linescanData = reshape(data,numChannels,samplesPerFrame,[]);
        linescanData = linescanData(1,:,:);
        linescanData = permute(linescanData,[2 3 1]);
        numLines = size(linescanData,2);
        channelAvg = [channelAvg  mean(linescanData,'omitnan')];
    end

    figure()
    plot(channelAvg)
    title(['m' num2str(m) ', c' num2str(c)])

    %save data
    save(fullfile(cur_dir,'Ch1_triggers_extracted.mat'),'channelAvg')
end



%% (ALREADY DONE) Pre-process raw velocity, fluorescence, width data

waitbarfig = waitbar(0,'Mouse 1, condition 1');
overwrite_vectors = false;
for m = 1:num_mice
    for c = 1:num_conds
        %preallocate space:
        %vessel arrays
        avg_vels = nan([num_vessels, num_pts]);
        responding_vels = nan([num_vessels, 1]);
        
        %neuron arrays
        avg_fluo = nan([num_neurons, num_pts]);
        avg_neuron = nan([num_neurons, 1]);
        responding_neurons_on = nan([num_neurons, 1]);
        responding_neurons_off = nan([num_neurons, 1]);
        
        %width arrays
        avg_widths = nan([num_vessels, num_pts]);
        responding_widths = nan([num_vessels, 1]);


        waitbar((m-1)/num_mice + (c-1)/(num_conds*num_mice),waitbarfig,['Mouse ' num2str(m) ', condition ' num2str(c) ': loading velocity data']);
        
        cur_dir = fullfile(base_folder,l(m).name,cond_names{c},'stim');
        if exist(cur_dir,'dir') && (overwrite_vectors || ~exist(fullfile(cur_dir,"vectors.mat"),'file'))

            %get vessel velocities
            for v = 1:num_vessels
                if foundvels(c,m,v)
                    curfilename = velnames{c,m,v};
    
                    %load velocity vector
                    load(curfilename)
                    
                    %load triggers
                    exp_folder = fileparts(curfilename);
                    curtrigfilename = ls(fullfile(exp_folder,'Ch1_triggers_extracted*'));
                    if ~isempty(curtrigfilename)
                        assert(size(curtrigfilename,1)==1,'unexpected # of trigger files')
                        curtrigfilename = fullfile(exp_folder,curtrigfilename);
                        load(curtrigfilename)
                        channelAvg = medfilt1(channelAvg,500);
                        minTrig = min(channelAvg);
                        maxTrig = max(channelAvg);
                        channelAvg = channelAvg>mean([minTrig maxTrig]);
                        risefall = find(diff(channelAvg)~=0);
                        trig_inds = round(risefall(1:2:end)/50);
                        num_reps = length(risefall)/2;
                        assert((num_reps==25)|(num_reps==20),'unexpected number of stimuli')
                        velocities = Result(:,3);
                        medvel = median(velocities);
                        if medvel<0
                            medvel = -medvel;
                            velocities = -velocities;
                        end
                        velocities(velocities>(1.6*medvel)) = nan;
                        velocities(velocities<(0.5*medvel)) = nan;
                        for i = 1:3
                            velocities_smooth = movmedian(velocities,20);
                            bad_idx = abs(velocities_smooth-velocities)>(0.2*velocities_smooth);
                            velocities(bad_idx) = nan;
                            velocities = fillmissing(velocities,"linear",'MissingLocations',isnan(velocities));
                        end
                        velocities_smooth = movmean(velocities,3);
    
                        %debugging:
                        % figure()
                        % plot(Result(:,3))
                        % hold on
                        % plot(velocities_smooth)
                        % title(['m' num2str(m) ' c' num2str(c) ' v' num2str(v)])
    
                        ONresponsevals = nan(num_reps,1);
                        tmps = nan(num_reps,length(time_vec));
                        for r = 1:num_reps
                            inds = (trig_inds(r)-10):(trig_inds(r)+num_pts-11);
                            inds(inds>length(velocities_smooth)) = [];
                            tmp = velocities_smooth(inds);
                            ONresponsevals(r) = mean(tmp(stim_inds)) - mean(tmp(baseline_inds));
                            tmps(r,1:length(tmp)) = tmp;
                        end
                        [~, h] = signrank(ONresponsevals);
                        if mean(ONresponsevals)<0 && h==1
                            h = -1;
                        end
                        responding_vels(v) = h;
                        avg_vels(v,:) = mean(tmps,'omitnan');
                    end
                end
            end
    
    
            %get vessel widths
            waitbar((m-1)/num_mice + (c-1)/(num_conds*num_mice),waitbarfig,['Mouse ' num2str(m) ', condition ' num2str(c) ': loading width data']);
            
            for v = 1:num_vessels
                if founddiams(c,m,v)
                    curfilename = diamnames{c,m,v};
    
                    %load diam vector
                    load(curfilename)
                    
                    %load triggers
                    exp_folder = fileparts(curfilename);
                    curtrigfilename = ls(fullfile(exp_folder,'Ch1_triggers_extracted*'));
                    if ~isempty(curtrigfilename)
                        assert(size(curtrigfilename,1)==1,'unexpected # of trigger files')
                        curtrigfilename = fullfile(exp_folder,curtrigfilename);
                        load(curtrigfilename)
                        channelAvg = medfilt1(channelAvg,500);
                        minTrig = min(channelAvg);
                        maxTrig = max(channelAvg);
                        channelAvg = channelAvg>mean([minTrig maxTrig]);
                        risefall = find(diff(channelAvg)~=0);
                        trig_inds = round(risefall(1:2:end)/50);
                        num_reps = length(risefall)/2;
                        assert((num_reps==25)|(num_reps==20),'unexpected number of stimuli')
                        widths = results.vesselDiameterInMicronsLineVec;
                        medwidth = median(widths);
                        if medwidth<0
                            medwidth = -medwidth;
                            widths = -widths;
                        end
                        widths(widths>(1.6*medwidth)) = nan;
                        widths(widths<(0.5*medwidth)) = nan;
                        for i = 1:3
                            widths_smooth = movmedian(widths,20);
                            bad_idx = abs(widths_smooth-widths)>(0.2*widths_smooth);
                            widths(bad_idx) = nan;
                            widths = fillmissing(widths,"linear",'MissingLocations',isnan(widths));
                        end
                        widths_smooth = movmean(widths,30);
    
                        %debugging widths
                        % figure()
                        % plot(results.vesselDiameterInMicronsLineVec)
                        % hold on
                        % plot(widths_smooth)
                        % title(['m' num2str(m) ' c' num2str(c) ' v' num2str(v)])
    
                        ONresponsevals = nan(num_reps,1);
                        tmps = nan(num_reps,length(time_vec));
                        for r = 1:num_reps
                            inds = (trig_inds(r)-10):(trig_inds(r)+num_pts-11);
                            inds(inds>length(widths_smooth)) = [];
                            tmp = widths_smooth(inds);
                            ONresponsevals(r) = mean(tmp(31:70),'omitnan') - mean(tmp(1:10),'omitnan');
                            tmps(r,1:length(inds)) = tmp;
                        end
                        [~, h] = signrank(ONresponsevals);
                        if mean(ONresponsevals)<0 && h==1
                            h = -1;
                        end
                        responding_widths(v) = h;
                        avg_widths(v,:) = mean(tmps,'omitnan');
                    end
                end
            end
    
    
            %get neural activity
            for n = 1:6
                if foundneurons(c,m,n)

                    waitbar((m-1)/num_mice + (c-1)/(num_conds*num_mice),waitbarfig,['Mouse ' num2str(m) ', condition ' num2str(c) ': loading neural activity data (neuron ' num2str(n) ')']);
                    
                    curfilename = neuronnames{c,m,n};
                    
                    %load tif file
                    cur_image = read_file(curfilename);
                    [h,w,nf] = size(cur_image);
    
                    %%collapse each linescan to a single fluorescence value (average of entire line)
                    cur_image(cur_image==0) = nan;
                    fluo_vec_1 = mean(cur_image,2,'omitnan');
                    fluo_vec_1 = reshape(fluo_vec_1,[1, numel(fluo_vec_1)]);
    
                    %shrink to match sample rate of velocities
                    fluo_vec_1 = imresize(fluo_vec_1,[1 length(velocities_smooth)],"bicubic"); 
                    
                    %get an average linescan, for later estimates of neuron cross section
                    cur_image = permute(cur_image,[2 1 3]); %[line, morelines, moreframes]
                    cur_image = reshape(cur_image,[w h*nf]); %[line, alllines]
                    fluo_vec_2 = median(cur_image,2,'omitnan'); %avg line
                    if size(avg_neuron,2)<w
                        avg_neuron(:,end+1:w) = nan;
                    end
                    avg_neuron(n,1:w) = fluo_vec_2;
    
                    %load triggers
                    if n==1
                        curtrigfilename = ls(fullfile(exp_folder,'Ch1_triggers_extracted*'));
                        if ~isempty(curtrigfilename)
                            assert(size(curtrigfilename,1)==1,'unexpected # of trigger files')
                            curtrigfilename = fullfile(exp_folder,curtrigfilename);
                            load(curtrigfilename)
                            channelAvg = medfilt1(channelAvg,500);
                            minTrig = min(channelAvg);
                            maxTrig = max(channelAvg);
                            channelAvg = channelAvg>mean([minTrig maxTrig]);
                            risefall = find(diff(channelAvg)~=0);
                            trig_inds = round(risefall(1:2:end)/50);
                            num_reps = length(risefall)/2;
                            assert((num_reps==25)|(num_reps==20),'unexpected number of stimuli')
                        end
                    end
                    if ~isempty(curtrigfilename)
                        %parse into stimulus timepoints 
                        ONresponsevals = nan(num_reps,1);
                        OFFresponsevals = nan(num_reps,1);
                        tmps = nan(num_reps,length(time_vec));
                        for r = 1:num_reps
                            inds = (trig_inds(r)-10):(trig_inds(r)+num_pts-11);
                            inds(inds>length(fluo_vec_1)) = [];
                            tmp = fluo_vec_1(inds);
                            ONresponsevals(r) = mean(tmp(stim_inds),'omitnan') - mean(tmp(baseline_inds),'omitnan');
                            OFFresponsevals(r) = mean(tmp(post_stim_inds),'omitnan') - mean(tmp(baseline_inds),'omitnan');
                            tmps(r,1:length(inds)) = tmp;
                        end
                        if ~all(isnan(ONresponsevals))
                            [~, h] = signrank(ONresponsevals);
                            if mean(ONresponsevals)<0 && h==1
                                h = -1;
                            end
                            responding_neurons_on(n) = h;
                        end
                        if ~all(isnan(OFFresponsevals))
                            [~, h] = signrank(OFFresponsevals);
                            if mean(OFFresponsevals)<0 && h==1
                                h = -1;
                            end
                            responding_neurons_off(n) = h;
                        end
                        avg_fluo(n,:) = mean(tmps,1,'omitnan');
                    end
                end
            end
    
            %end mouse/condition loop
            save(fullfile(cur_dir,"vectors.mat"),'avg_neuron','avg_fluo','responding_neurons_on','responding_neurons_off','avg_vels','responding_vels','avg_widths','responding_widths');
    
        end
    end
end
delete(waitbarfig)


%% (STEP 2) load all processed data

% load all vectors into a single array
%vessel arrays
avg_vels_all = nan([num_conds, num_mice, num_vessels, num_pts]);
responding_vels_all = nan([num_conds, num_mice, num_vessels]);

%neuron arrays
avg_fluo_all = nan([num_conds, num_mice, num_neurons, num_pts]);
avg_neuron_all = nan([num_conds, num_mice, num_neurons, 1]);
responding_neurons_on_all = nan([num_conds, num_mice, num_neurons]);
responding_neurons_off_all = nan([num_conds, num_mice, num_neurons]);

%width arrays
avg_widths_all = nan([num_conds, num_mice, num_vessels, num_pts]);
responding_widths_all = nan([num_conds, num_mice, num_vessels]);

%loop through directories to load processed data vectors
for m = 1:num_mice
    for c = 1:num_conds
        cur_dir = fullfile(base_folder,l(m).name,cond_names{c},'stim');
        if exist(fullfile(cur_dir,"vectors.mat"),'file')
            
            load(fullfile(cur_dir,"vectors.mat"))

            avg_vels_all(c,m,:,:) = avg_vels;
            responding_vels_all(c,m,:) = responding_vels;

            avg_fluo_all(c,m,:,:) = avg_fluo;
            avg_neuron = squeeze(avg_neuron);
            ln = size(avg_neuron,2);
            if size(avg_neuron_all,4)<ln
                avg_neuron_all(:,:,:,end+1:ln) = nan;
            end
            avg_neuron_all(c,m,:,1:ln) = avg_neuron;
            responding_neurons_on_all(c,m,:) = responding_neurons_on;
            responding_neurons_off_all(c,m,:) = responding_neurons_off;
            
            avg_widths_all(c,m,:,:) = avg_widths;
            responding_widths_all(c,m,:) = responding_widths;

        end
    end
end

%shorten names:
avg_vels = avg_vels_all;
responding_vels = responding_vels_all;
avg_fluo = avg_fluo_all;
avg_neuron = avg_neuron_all;
responding_neurons_on = responding_neurons_on_all;
responding_neurons_off = responding_neurons_off_all;
avg_widths = avg_widths_all;
responding_widths = responding_widths_all;


%organize into long data format
%vessel vectors: [vessel, datatype, timepoint]
%vessel values: [vessel, datatype]
%vessel datatypes
% 1: velocity
% 2: dV/V
% 3: respondingvel
% 4: width
% 5: dW/W
% 6: respondingwidth (-1/0/1)
% 7: mouseID (1-27)
% 8: vesselID (1-3)
% 9: sex (1=male, 2=female)
% 10: condition (8 conds, but 4 for deltas: 1=control, 2=psilo, 3=mdl, 4=psilo+mdl)

%neuron vectors: [neuron, datatype, timepoint]
%neuron values: [neuron, datatype]
%neuron datatypes
% 1: F
% 2: dF/F
% 3: respondingon
% 4: respondingoff
% 5: mouseID (1-27)
% 6: neuronID (1-6)
% 7: sex (1=male, 2=female)
% 8: condition (8 conds, but 4 for deltas: 1=control, 2=psilo, 3=mdl, 4=psilo+mdl)

max_vessels_long = num_conds*num_mice*num_vessels;
max_neurons_long = num_conds*num_mice*num_neurons;
vessel_long_vector = nan(max_vessels_long,10,num_pts);
neuron_long_vector = nan(max_neurons_long,8,num_pts);
cur_vessel = 0;
cur_neuron = 0;
%widefield:
% sex_dictionary = dictionary(...
%     {'1L','2R','3N','4B','3108','3113','3114','3115','3815','3817','3904','3905','6186','6187','6189','543959','546842','547880'},...
%     [ 1    1    1    1    1      2      2      2      2      2      2      2      2      2      2      2        1        1      ]);
%2P:
sex_dictionary = dictionary(...
    {'1L','2R','3N','4B','885','2114','2117','2119','2236','2621','3108','3109','3113','3114','3115','3815','3817','3904','3905','6184N_M','6185R_M','6186','6187_F','6189','543959','546842','547880', '6680_None', '6681_Top', '6682_Bottom', '7675', '7676'},...
    [ 1    1    1    1    2     1      1      2      1      1      1      1      2      2      2      2      2       2     2      1         1         2      2        2      2        1        1         1            1           1              2       2   ]);

%loop through data to organize into long data format
for m = 1:num_mice
    for v = 1:num_vessels
        for c = 1:num_conds
            cur_vessel = cur_vessel + 1;
            vessel_long_vector(cur_vessel,1,:) = avg_vels(c,m,v,:);
            vessel_long_vector(cur_vessel,2,:) = 100*((avg_vels(c,m,v,:)/mean(avg_vels(c,m,v,baseline_inds),'omitnan'))-1);
            vessel_long_vector(cur_vessel,3,:) = responding_vels(c,m,v);
            vessel_long_vector(cur_vessel,4,:) = avg_widths(c,m,v,:);
            vessel_long_vector(cur_vessel,5,:) = 100*((avg_widths(c,m,v,:)/mean(avg_widths(c,m,v,baseline_inds),'omitnan'))-1);
            vessel_long_vector(cur_vessel,6,:) = responding_widths(c,m,v);
            vessel_long_vector(cur_vessel,7,:) = m;
            vessel_long_vector(cur_vessel,8,:) = v;
            vessel_long_vector(cur_vessel,9,:) = sex_dictionary({l(m).name});
            vessel_long_vector(cur_vessel,10,:) = c;
        end
    end

    for n = 1:num_neurons
        for c = 1:num_conds
            cur_neuron = cur_neuron + 1;
            neuron_long_vector(cur_neuron,1,:) = avg_fluo(c,m,n,:);
            neuron_long_vector(cur_neuron,2,:) = 100*((avg_fluo(c,m,n,:)/mean(avg_fluo(c,m,n,baseline_inds),'omitnan'))-1);
            neuron_long_vector(cur_neuron,3,:) = responding_neurons_on(c,m,n);
            neuron_long_vector(cur_neuron,4,:) = responding_neurons_off(c,m,n);
            neuron_long_vector(cur_neuron,5,:) = m;
            neuron_long_vector(cur_neuron,6,:) = n;
            neuron_long_vector(cur_neuron,7,:) = sex_dictionary({l(m).name});
            neuron_long_vector(cur_neuron,8,:) = c;
        end
    end
end

save("C:\Users\misaa\Desktop\neuron_long_vector.mat","neuron_long_vector")
save("C:\Users\misaa\Desktop\vessel_long_vector.mat","vessel_long_vector")



%% (STEP 3) exclude bad data (found from next section), create inclusive categories
%neurons and vessel data are only excluded if by manual inspection it is
%noted that the data quality is very poor. For velocity, sometimes the
%radon transform does not confidently find velocity measurements leading to
%large missing data caps (too large to be reasonably interpolated). For
%neurons, sometimes baseline fluorescence is too noisy or too small,
%leading to undreliable dF/F. For width, sometimes the Gaussian function
%(or two term function) does not fit well to the cross-section profile, or
%the sampling is clipped too much on one edge of the vessel, or the fitted
%guassian is too large to reasonably be a capillary. Details on the reason
%for each excluded datapoint is in the next section.

%exlude format: [mouse cond vessel/neuron]
bad_vels_inds = [5 1 1;...
    7 3 3;...
    13 1 1;...
    14 2 1;...
    24 7 2]; 

bad_widths_inds = [6 3 2;...
    14 1 3;...
    14 2 3;...
    14 3 3;...
    17 1 3;...
    17 2 3;...
    14 4 3;...
    21 2 3;...
    21 3 3;...
    21 4 3;...
    21 7 3;...
    21 8 3;...
    23 5 3;...
    23 6 3;...
    24 7 2;...
    24 8 2];

bad_neurons_inds = [7 1 1;...
    7 2 1;...
    13 1 5;...
    13 2 5;...
    13 3 5;...
    13 3 1;...
    13 3 2;...
    13 3 3;...
    13 3 4;...
    13 3 5;...
    13 3 6;...
    13 4 1;...
    13 4 2;...
    13 4 3;...
    13 4 4;...
    13 4 5;...
    13 4 6];

c2s = [2 1 4 3 6 5 8 7]; %complementary conditions for c=[1 2 3 4 5 6 7 8]

% exclude bad velocities
for i = 1:size(bad_vels_inds,1)
    m = bad_vels_inds(i,1);
    v = bad_vels_inds(i,3);
    c1 = bad_vels_inds(i,2);
    c2 = c2s(c1);

    %delete bad velocity data
    rows_idx = find(vessel_long_vector(:,7)==m & vessel_long_vector(:,8)==v & (vessel_long_vector(:,10)==c1 | vessel_long_vector(:,10)==c2));
    vessel_long_vector(rows_idx,1:3,:) = nan;
end

% exclude bad widths
for i = 1:size(bad_widths_inds,1)
    m = bad_widths_inds(i,1);
    v = bad_widths_inds(i,3);
    c1 = bad_widths_inds(i,2);
    c2 = c2s(c1);

    %delete bad velocity data from long vectors
    rows_idx = find(vessel_long_vector(:,7)==m & vessel_long_vector(:,8)==v & (vessel_long_vector(:,10)==c1 | vessel_long_vector(:,10)==c2));
    vessel_long_vector(rows_idx,4:6,:) = nan;
end


% exclude bad neurons
for i = 1:size(bad_neurons_inds,1)
    m = bad_neurons_inds(i,1);
    n = bad_neurons_inds(i,3);
    c1 = bad_neurons_inds(i,2);
    c2 = c2s(c1);

    %delete bad velocity data from long vectors
    rows_idx = find(neuron_long_vector(:,5)==m & neuron_long_vector(:,6)==n & (neuron_long_vector(:,8)==c1 | neuron_long_vector(:,8)==c2));
    neuron_long_vector(rows_idx,1:4,:) = nan;
end


% create vessel and neuron categories
max_neurons = size(neuron_long_vector,1);
max_vessels = size(vessel_long_vector,1);

%list of categories
%all = all neurons/vessels that were measured in both pre and post conditions
%pos = vessels that had a positive response during the stimulus
%neg = vessels that had a negative response during the stimulus
%none = vessels that had no significant response during the stimulus
%on = neurons positively responsive during stimulus (but not during offset) in each timepoint
%off = neurons positively responsive during offset (but not during stimulus) in each timepoint
%mixed = neurons positively responsive during stimulus and offset each timepoint
%none = neurons not responsive in each timepoint
%either_"" = neurons/vessels that met the condition in either the pre or post (or both) timepoint

% vessel_long_vector: [vessel, datatype, timepoint]
%vessel datatypes
% 1: velocity
% 2: dV/V
% 3: respondingvel
% 4: width
% 5: dW/W
% 6: respondingwidth (-1/0/1)
% 7: mouseID (1-27)
% 8: vesselID (1-3)
% 9: sex (1=male, 2=female)
% 10: condition (8 conds)

% neuron_long_vector: [neuron, datatype, timepoint]
%neuron datatypes
% 1: F (raw)
% 2: dF/F
% 3: respondingon
% 4: respondingoff
% 5: mouseID (1-27)
% 6: neuronID (1-6)
% 7: sex (1=male, 2=female)
% 8: condition (8 conds)


%classify/categorize neurons and vessels
all_vels = ~isnan(vessel_long_vector(:,3,1)); %every row that's not a nan
all_vels(1:2:end) = all_vels(1:2:end) & all_vels(2:2:end); %every non-nan pre row that has a non-nan post row
all_vels(2:2:end) = all_vels(1:2:end); %every non-nan post row that has a non-nan pre row
pos_vels = all_vels & vessel_long_vector(:,3,1)==1;
neg_vels = all_vels & vessel_long_vector(:,3,1)==-1;
none_vels = all_vels & vessel_long_vector(:,3,1)==0;

all_widths = ~isnan(vessel_long_vector(:,6,1)); %every row that's not a nan
all_widths(1:2:end) = all_widths(1:2:end) & all_widths(2:2:end); %every non-nan pre row that has a non-nan post row
all_widths(2:2:end) = all_widths(1:2:end); %every non-nan post row that has a non-nan pre row
pos_widths = all_vels & vessel_long_vector(:,6,1)==1;
neg_widths = all_vels & vessel_long_vector(:,6,1)==-1;
none_widths = all_vels & vessel_long_vector(:,6,1)==0;


all_neurons = ~isnan(neuron_long_vector(:,3,1)); %every row that's not a nan
all_neurons(1:2:end) = all_neurons(1:2:end) & all_neurons(2:2:end); %every non-nan pre row that has a non-nan post row
all_neurons(2:2:end) = all_neurons(1:2:end); %every non-nan post row that has a non-nan pre row
pos_neurons = all_neurons & neuron_long_vector(:,3,1)==1;
neg_neurons = all_neurons & neuron_long_vector(:,3,1)==-1;
none_neurons = all_neurons & neuron_long_vector(:,3,1)==0;
%classify neurons based on their classifications in either pre or post
either_pos_neurons = pos_neurons; 
either_pos_neurons(1:2:end) = either_pos_neurons(1:2:end) | either_pos_neurons(2:2:end);
either_pos_neurons(2:2:end) = either_pos_neurons(1:2:end);
either_neg_neurons = neg_neurons; 
either_neg_neurons(1:2:end) = either_neg_neurons(1:2:end) | either_neg_neurons(2:2:end);
either_neg_neurons(2:2:end) = either_neg_neurons(1:2:end);
either_none_neurons = none_neurons; 
either_none_neurons(1:2:end) = either_none_neurons(1:2:end) | either_none_neurons(2:2:end);
either_none_neurons(2:2:end) = either_none_neurons(1:2:end);
%classify neurons by their pre-treatment classification
pre_pos_neurons = pos_neurons; 
pre_pos_neurons(2:2:end) = pre_pos_neurons(1:2:end);
pre_neg_neurons = neg_neurons; 
pre_neg_neurons(2:2:end) = pre_neg_neurons(1:2:end);
pre_none_neurons = none_neurons; 
pre_none_neurons(2:2:end) = pre_none_neurons(1:2:end);


%% (OPTIONAL) inspect data
%create 3 plots for all data of the specified mice, conditions, and neurons/vessels (dF/F, dW/W, dV/V)
m = 1;
c = 1:8;
v = 1:3;
n = 1:6;

figure()

%subplot 1: inspect all neuron traces
rows_to_plot = find(all_neurons & any(neuron_long_vector(:,5,1)==m,2) & any(neuron_long_vector(:,6,1)==n,2) & any(neuron_long_vector(:,8,1)==c,2));
subplot(3, 1, 1)
tmp = squeeze(neuron_long_vector(rows_to_plot,2,:))';
tmp(isnan(tmp)) = 0;
plot(time_vec,tmp)
title('dF/F')

%subplot 2: inspect all widths
rows_to_plot = find(all_widths & any(vessel_long_vector(:,7,1)==m,2) & any(vessel_long_vector(:,8,1)==v,2) & any(vessel_long_vector(:,10,1)==c,2));
subplot(3, 1, 2)
tmp = squeeze(vessel_long_vector(rows_to_plot,5,:))';
tmp(isnan(tmp)) = 0;
plot(time_vec,tmp)
title('dW/W')

%subplot 3: inspect all vels
rows_to_plot = find(all_vels & any(vessel_long_vector(:,7,1)==m,2) & any(vessel_long_vector(:,8,1)==v,2) & any(vessel_long_vector(:,10,1)==c,2));
subplot(3, 1, 3)
tmp = squeeze(vessel_long_vector(rows_to_plot,2,:))';
tmp(isnan(tmp)) = 0;
plot(time_vec,tmp)
title('dV/V')

%bad data:
%mouse 5: c1v1 bad vel - missing datapoints
%mouse 6: c3v2 bad width - bad fit
%mouse 7: c3v3 bad vel; c1/2,n1 missing neurons
%mouse 13: c1,v1 bad vel data; c1-3,n5 missing neurons
%mouse 13: c3/4: neuron baselines approaching 0 (not good for dF/F, exclude)
%mouse 14: c2,v1 missing vel; c1-3,v3 missing width
%mouse 17: c1/2/4,v3 no good width measurement
%mouse 21: c2/3/4/7/8,v3 no good width measurement
%mouse 23: c5/6,v3 no good width measurement
%mouse 24: c7,v2 missing vel; c7-8,v2 no good width measurement


%% (OPTIONAL) check distribution of neuron and vessel baselines (so that dX/X is reasonable)

figure()
subplot(1,3,1)
baselines = mean(neuron_long_vector(:,1,baseline_inds),3,'omitnan');
baselines = baselines(all_neurons);
histogram(baselines)
title(['F (min ' num2str(min(baselines)) ')'])

subplot(1,3,2)
baselines = mean(vessel_long_vector(:,1,baseline_inds),3,'omitnan');
baselines = baselines(all_vels);
histogram(baselines)
title(['V (min ' num2str(min(baselines)) ')'])

subplot(1,3,3)
baselines = mean(vessel_long_vector(:,4,baseline_inds),3,'omitnan');
baselines = baselines(all_widths);
histogram(baselines)
title(['W (min ' num2str(min(baselines)) ')'])


%% count number of samples/mice per dataset
counts = nan([4 9]);
counts(:,1) = [1 3 5 7];
for c = 1:4
    cond = counts(c,1);

    %on neurons
    rows = find(either_pos_neurons & neuron_long_vector(:,8,1)==cond);
    counts(c,2) = length(rows); %count vessels in this condition
    mice = neuron_long_vector(rows,5,1);
    counts(c,3) = length(unique(mice)); %count mice in this condition

    %all neurons
    rows = find(all_neurons & neuron_long_vector(:,8,1)==cond);
    counts(c,4) = length(rows); %count vessels in this condition
    mice = neuron_long_vector(rows,5,1);
    counts(c,5) = length(unique(mice)); %count mice in this condition

    %vels
    rows = find(all_vels & vessel_long_vector(:,10,1)==cond);
    counts(c,6) = length(rows); %count vessels in this condition
    mice = vessel_long_vector(rows,7,1);
    counts(c,7) = length(unique(mice)); %count mice in this condition

    %widths
    rows = find(all_widths & vessel_long_vector(:,10,1)==cond);
    counts(c,8) = length(rows); %count vessels in this condition
    mice = vessel_long_vector(rows,7,1);
    counts(c,9) = length(unique(mice)); %count mice in this condition
end
count_table = array2table(counts,'VariableNames',{'condition','pos neurons','pos neuron-mice','all neurons','all neuron-mice','velocity-vessels','velocity-mice','width-vessels','width-mice'})



%% set basic figure panel parameters
%list of anel types
%raster plot (6 total: control/psilo/psilo+mdl, neuron/velocity)
%mean dx/x +/- sd/sem (3 total: neuronON/velocity/width)
%delta dx/x +/- sd/sem (3 total: neuronON/velocity/width)
%distribution of responsivitiy types (3 total: neuron/velocity/width)
%boxplots of baselines (9 total: control/psilo/psilo+mdl, neuronON/velocity/width)
%boxplots of deltas (9 total: control/psilo/psilo+mdl, neuronON/velocity/width)

%set general figure parameters

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
ylabels_baselines = {'raw brightness','velocity (\mm/s)','width (um)'};
font_size = 4.8; %[4 4.8 5.6]Matlab = [5 6 7]Inkscape

vector_plot_size = [150 85];
box_plot_size = [120 85];
raster_plot_size = [130 250];


%% (Plotting) create raster plots
plotlabels_c2 = {'pre','post'};
plotlabels_c = {'control','','psilocybin','','MDL','','psilocybin+MDL'};
ylabels_dt = {'neuron','vessel','vessel'};
raster_ylimits = [-40 60; -20 20; -3 3];

for dt = [1 2] %neuron, velocity
    for c = [7] %control pre/post, psilo pre/post, psilo+mdl pre/post
        figure('Position',[100 100 raster_plot_size])

        %get sorted order of rows rows
        switch dt
            case 1
                rows = find(all_neurons & neuron_long_vector(:,8,1)==c); %get all neurons for this condition
                %sort by increasing mean responses (stim - post-stim)
                sort_values = mean(neuron_long_vector(rows,2,stim_inds),3,'omitnan') - mean(neuron_long_vector(rows,2,post_stim_inds),3,'omitnan');
                [~, order] = sort(sort_values);
                rows = rows(order);
            case 2 
                rows = find(all_vels & vessel_long_vector(:,10,1)==c); %get all vessels for this condition
                %sort by increasing mean response (during stim)
                sort_values = mean(vessel_long_vector(rows,2,stim_inds),3,'omitnan');
                [~, order] = sort(sort_values);
                rows = rows(order);
            case 3
                rows = find(all_widths & vessel_long_vector(:,10,1)==c); %get all vessels for this condition
                %sort by increasing mean response (during stim)
                sort_values = mean(vessel_long_vector(rows,5,stim_inds),3,'omitnan');
                [~, order] = sort(sort_values);
                rows = rows(order);
        end

        for c2 = 1:2 %loop for pre/post condition
            %get data for current condition
            switch dt
                case 1
                    sorted_vectors = squeeze(neuron_long_vector(rows-1+c2,2,:));
                case 2 
                    sorted_vectors = squeeze(vessel_long_vector(rows-1+c2,2,:));
                case 3
                    sorted_vectors = squeeze(vessel_long_vector(rows-1+c2,5,:));
            end

            %create raster plot for this datatype/condition
            axes('InnerPosition',[0.12+(c2-1)*0.4 0.1 0.35 0.7])
            % plot(sorted_vectors')
            for row = 1:length(rows)
                for s = 1:(num_pts-1)
                    val = sorted_vectors(row,s);
                    if val>=0
                        colorval = min([1, val/raster_ylimits(dt,2)]);
                        cur_color = (SN_colors(5,:)*colorval + brightwhite*(1-colorval)); %0=white, 1=pure color
                    else
                        colorval = min([1, val/raster_ylimits(dt,1)]);
                        cur_color = (SN_colors(10,:)*colorval + brightwhite*(1-colorval)); %0=white, 1=pure color
                    end
                    patch(time_vec([s s+1 s+1 s]),[row-0.5 row-0.5 row+0.5 row+0.5],'k','EdgeColor','none','FaceColor',cur_color)
                end
            end

            %plot stimulus region
            patch([-4 20 20 -4],[-0.01 -0.01 -0.05 -0.05]*length(rows),'w','EdgeColor','none')
            patch([0 6 6 0],[-0.01 -0.01 -0.04 -0.04]*length(rows),'k','EdgeColor','none','FaceColor',[0.3 0.3 0.3])
            text(3,-0.025*length(rows),'stim','FontSize',font_size,'FontName','Arial','HorizontalAlignment','center','VerticalAlignment','middle','Color',brightwhite)
            
            %add subplot elements
            ax = gca;
            ax.YDir = 'reverse';
            ax.YLim = [-0.05*length(rows) length(rows)+0.5];
            ax.XLim = [-3 15];
            text(6,-0.09*length(rows),[plotlabels_c2{c2} ' ' plotlabels_c{c}],'FontSize',font_size,'FontName','Arial','HorizontalAlignment','center','Color',[0 0 0])
            ax.FontSize = font_size;
            ax.FontName = 'Arial';
            if c2==1
                ax.YTick = [1 length(rows)];
                ylabel(ylabels_dt{dt},'FontSize',font_size)
                xlabel('time (s)','FontSize',font_size)
            else
                ylabel('')
                ax.YTick = [];
            end
        end

        print(gcf,'-vector','-dsvg',fullfile(base_folder,['raster ' ylabels_dt{dt} ' ' num2str(c) '.svg']))

    end
end


%% (Plotting) plot distribution of baseline fluorescence, flow speeds, widths
box_ylimits = [0 0; -1 2; -20 30];
cur_comparison = 0;
lme_baseline_ps = nan(2,1);
dtypenames = {'','velocities','widths'};
xlabels_histograms = {'','velocity (mm/s)' 'width (um)'};
ylabels_histograms = {'','number of vessels', 'number of vessels'};
xlimits_histograms = [0 0; 0 2; 0 30];
ylimits_histograms = [0 0; 0 25; 0 25];
numbins = [0 20 30];
% conds_to_compare = [2 4];  %control pre/post, psilo pre/post
conds_to_compare = [4 8];  %control pre/post, psilo pre/post
boxplot_colors = [colors(conds_to_compare(1)-1,:); colors(conds_to_compare(1),:); colors(conds_to_compare(1),:); colors(conds_to_compare(2)-1,:); colors(conds_to_compare(2),:); colors(conds_to_compare(2),:)];
for dt = [2 3] %velocity, width (skip fluorescence - not meaningful)
    %allocate space for baseline values
    baseline_vals = nan(1,8);
    Data = [];
    Treatment = [];
    Mouse = [];
    Timepoint = [];
    Cell = [];
    % for c = [conds_to_compare(1)-1 conds_to_compare(1) conds_to_compare(2)-1 conds_to_compare(2) conds_to_compare(3)-1 conds_to_compare(3)]
    for c = [conds_to_compare(1)-1 conds_to_compare(1) conds_to_compare(2)-1 conds_to_compare(2)]
        switch dt
            case 2 
                rows = find(all_vels & vessel_long_vector(:,10,1)==c); %get all vessels for this condition
                nr = length(rows);
                cur_cells = vessel_long_vector(rows,8,1);
                mean_vec = squeeze(mean(vessel_long_vector(rows,1,:),1,'omitnan'))'; %mean V vector
                std_vec = squeeze(std(vessel_long_vector(rows,1,:),0,1,'omitnan'))';
                sem_vec = std_vec/sqrt(nr);
                cur_baseline_vals = squeeze(mean(vessel_long_vector(rows,1,baseline_inds),3,'omitnan')); %mean V during baseline
                cur_mouse_vals = vessel_long_vector(rows,7); %store for long format
            case 3
                rows = find(all_widths & vessel_long_vector(:,10,1)==c); %get all vessels for this condition
                nr = length(rows);
                cur_cells = vessel_long_vector(rows,8,1);
                mean_vec = squeeze(mean(vessel_long_vector(rows,4,:),1,'omitnan'))'; %mean W vector
                std_vec = squeeze(std(vessel_long_vector(rows,4,:),0,1,'omitnan'))';
                sem_vec = std_vec/sqrt(nr);
                cur_baseline_vals = squeeze(mean(vessel_long_vector(rows,4,baseline_inds),3,'omitnan')); %mean W during baseline
                cur_mouse_vals = vessel_long_vector(rows,7); %store for long format
        end

        %store baselines
        if nr>size(baseline_vals,1)
            baseline_vals(end+1:nr,:) = nan;
        end
        baseline_vals(1:nr,c) = cur_baseline_vals;

        %store delta baseline and metadata in long format for anova
        if mod(c,2)==1
            pre_baseline_vals = cur_baseline_vals;
        else
            %store measurement and metadata in long format for stats
            Data = [Data; pre_baseline_vals];
            Mouse = [Mouse; cur_mouse_vals];
            Treatment = [Treatment; c*ones(nr,1)];
            Timepoint = [Timepoint; ones(nr,1)];
            Cell = [Cell; cur_cells];
    
            Data = [Data; cur_baseline_vals];
            Mouse = [Mouse; cur_mouse_vals];
            Treatment = [Treatment; c*ones(nr,1)];
            Timepoint = [Timepoint; 2*ones(nr,1)];
            Cell = [Cell; cur_cells];
        end
    end

    %create figure for baseline boxplots
    delta_baseline_vals = baseline_vals(:,conds_to_compare) - baseline_vals(:,conds_to_compare-1);
    % all_baseline_vals = [baseline_vals(:,[conds_to_compare(1)-1 conds_to_compare(1)]) delta_baseline_vals(:,1) nan(size(delta_baseline_vals,1),1) baseline_vals(:,[conds_to_compare(2)-1 conds_to_compare(2)]) delta_baseline_vals(:,2)];
    all_baseline_vals = [baseline_vals(:,[conds_to_compare(1)-1 conds_to_compare(1)]) delta_baseline_vals(:,1) nan(size(delta_baseline_vals,1),1) baseline_vals(:,[conds_to_compare(2)-1 conds_to_compare(2)]) delta_baseline_vals(:,2)];
    
    figure('Position',[100 100 box_plot_size])
    plot([0 8],[0 0],'--k')
    hold on
    x = kron([1 2 3 5 6 7],ones(size(all_baseline_vals,1),1));
    for j = 1:6
        swarmchart(x(:,j),all_baseline_vals(:,x(1,j)),1,0.8*boxplot_colors(j,:),'filled','o','MarkerFaceAlpha',1,'MarkerEdgeColor','none')
    end
    bp = boxplot(all_baseline_vals,'Symbol','');
    set(bp(6,:), 'Color', 'k'); %median
    set(bp(5,:), 'Color', 'k'); %box
    patch([3.75 4.25 4.25 3.75],[box_ylimits(dt,1)*1.1 box_ylimits(dt,:)*1.1 box_ylimits(dt,2)*1.1],'w','EdgeColor','None')
    plot([3 7],1.05*[box_ylimits(dt,2) box_ylimits(dt,2)],'k')
    plot([3 3],[0.95*box_ylimits(dt,2) 1.05*box_ylimits(dt,2)],'k')
    plot([7 7],[0.95*box_ylimits(dt,2) 1.05*box_ylimits(dt,2)],'k')
    text(3,1.1*box_ylimits(dt,1),'\Delta','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','center')
    text(7,1.1*box_ylimits(dt,1),'\Delta','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','center')
    ylim(box_ylimits(dt,:))
    xlim([0 8])
    set(gca,'FontSize',font_size)
    set(gca,'FontName','Arial')
    xticklabels({'pre','post','','','pre','post',''})
    ylim(box_ylimits(dt,:))
    ylabel(ylabels_baselines{dt},'FontSize',font_size)
    ax = gca;
    ax.Clipping = "off";
    box off
    
    %run LME for this dataset
    cur_comparison = cur_comparison + 1;
    lmedata.measurement = Data;
    lmedata.mouse = categorical(Mouse);
    lmedata.treatment = categorical(Treatment);
    lmedata.timepoint = categorical(Timepoint);
    lmedata.cell = categorical(Cell);
    lmetable = struct2table(lmedata);
    % treatment = fixed effect, mouse = random intercept
    % lme_model = fitlme(lmetable, 'measurement ~ treatment + (1|mouse)');
    lme_model = fitlme(lmetable, 'measurement ~ treatment*timepoint + (1|mouse) + (1|cell:mouse)');
    % lme_model = fitlme(lmetable, 'measurement ~ treatment*timepoint + order + (1|mouse) + (1|session:mouse) + (1|cell:session:mouse)');
            
    % ANOVA p-value for treatment effect
    % anova_results = anova(lme_model); 
    lme_baseline_ps(cur_comparison,:) = lme_model.Coefficients.pValue(4);


    %put placeholder significance bar
    text(5,1.1*box_ylimits(dt,2),'n.s.','FontSize',font_size,'VerticalAlignment','bottom','HorizontalAlignment','center')
    % text(5.7,1.09*box_ylimits(dt,2),['(p=' num2str(p(1),'%.02f') ')'],'FontSize',font_size,'VerticalAlignment','bottom','HorizontalAlignment','left')
    text(5.7,1.09*box_ylimits(dt,2),'(p=?','FontSize',font_size,'VerticalAlignment','bottom','HorizontalAlignment','left')

    print(gcf,'-vector','-dsvg',['C:\Users\misaa\Desktop\2p baseline boxplot conds ' num2str(conds_to_compare) ' dt ' dtypenames{dt} '.svg'])


    %create histsogram of all baseline values
    figure('Position',[100 100 box_plot_size])
    tmp = baseline_vals(:,conds_to_compare-1);
    plot([0 8],[0 0],'--k')
    hold on
    histogram(tmp(:),linspace(xlimits_histograms(dt,1),xlimits_histograms(dt,2),numbins(dt)),'FaceColor',colors(1,:))
    xlim(xlimits_histograms(dt,:))
    ylim(ylimits_histograms(dt,:))
    set(gca,'FontSize',font_size)
    set(gca,'FontName','Arial')
    xlabel(xlabels_histograms{dt},'FontSize',font_size)
    ylabel(ylabels_histograms{dt},'FontSize',font_size)
    box off

    print(gcf,'-vector','-dsvg',['C:\Users\misaa\Desktop\2p baseline hist conds ' num2str(conds_to_compare) ' dt ' dtypenames{dt} '.svg'])

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
    title('Residuals vs. Fitted');
    subplot(2,2,4);
    normplot(lme_residuals);
    title('normal probability plot');
end
lme_baseline_ps


%% (OPTIONAL) find datapoints with large baseline vessel widths
rows = find(Data>10);
big.baseline_width = Data(rows);
big.mouse_name = mouse_names(Mouse(rows));
big.treatment_name = treatment_names(Treatment(rows))';
big.timepoint_name = timepoint_names(Timepoint(rows))';
big.vessel_num = width_roi_nums(Cell(rows))';
big_width_table = struct2table(big);



%% (Plotting) create mean delta vectors
ylabels = {'pos neuron % \DeltaF/F','neg neuron % \DeltaF/F','all vessel % \DeltaV/V','all vessel % \DeltaW/W'};
vector_ylimits = [-5 15; -15 20; -10 10; -2 2]; %[neuron; velocity; width]
dtypenames = {'posneurons','negneurons','velocities','widths'};
% conds_to_compare = [1 2 3 4];
conds_to_compare = [7 8];
%for BOLD modeling:
vel_means = nan(length(time_vec),4); %[pts, cond(4)]

% vector_plot_size = [100 85];

for dt = [1 2 3 4] %posneuron, negneuron, velocity, width
    %allocate space for baseline values
    baseline_vals = nan(1,8);
    Data = [];
    Treatment = [];
    Mouse = [];

    %create figure for vectors
    % if dt<4
        % figure('Position',[100 100 vector_plot_size])
    % else
        figure('Position',[100 100 box_plot_size]) %for ED fig, shrink size
    % end
    
    plot([0 20],[0 0],'--k')
    hold on
    patch([0 6 6 0],[-100 -100 100 100],'k','EdgeColor','none','FaceAlpha',0.1)
    text(3,vector_ylimits(dt,2),'stimulus','Color',[0.2 0.2 0.2],'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',font_size)
    for c = conds_to_compare
        switch dt
            case 1
                rows = find(pre_pos_neurons & neuron_long_vector(:,8,1)==c); %get all (either_pos) neurons for this condition
                nr = length(rows);
                mean_vec = squeeze(mean(neuron_long_vector(rows,2,:),1,'omitnan'))'; %mean dF/F vector
                std_vec = squeeze(std(neuron_long_vector(rows,2,:),0,1,'omitnan'))';
                sem_vec = std_vec/sqrt(nr);
                cur_baseline_vals = squeeze(mean(neuron_long_vector(rows,1,baseline_inds),3,'omitnan')); %mean F during baseline
                cur_mouse_vals = neuron_long_vector(rows,5); %store for long format
            
            case 2
                rows = find(pre_neg_neurons & neuron_long_vector(:,8,1)==c); %get all (either_neg) neurons for this condition
                nr = length(rows);
                mean_vec = squeeze(mean(neuron_long_vector(rows,2,:),1,'omitnan'))'; %mean dF/F vector
                std_vec = squeeze(std(neuron_long_vector(rows,2,:),0,1,'omitnan'))';
                sem_vec = std_vec/sqrt(nr);
                cur_baseline_vals = squeeze(mean(neuron_long_vector(rows,1,baseline_inds),3,'omitnan')); %mean F during baseline
                cur_mouse_vals = neuron_long_vector(rows,5); %store for long format
            
            case 3
                rows = find(all_vels & vessel_long_vector(:,10,1)==c); %get all vessels for this condition
                nr = length(rows);
                mean_vec = squeeze(mean(vessel_long_vector(rows,2,:),1,'omitnan'))'; %mean dV/V vector
                std_vec = squeeze(std(vessel_long_vector(rows,2,:),0,1,'omitnan'))';
                sem_vec = std_vec/sqrt(nr);
                cur_baseline_vals = squeeze(mean(vessel_long_vector(rows,1,baseline_inds),3,'omitnan')); %mean V during baseline
                cur_mouse_vals = vessel_long_vector(rows,7); %store for long format

                vel_means(:,c) = mean_vec;
            case 4
                rows = find(all_widths & vessel_long_vector(:,10,1)==c); %get all vessels for this condition
                nr = length(rows);
                mean_vec = squeeze(mean(vessel_long_vector(rows,5,:),1,'omitnan'))'; %mean dW/W vector
                std_vec = squeeze(std(vessel_long_vector(rows,5,:),0,1,'omitnan'))';
                sem_vec = std_vec/sqrt(nr);
                cur_baseline_vals = squeeze(mean(vessel_long_vector(rows,4,baseline_inds),3,'omitnan')); %mean W during baseline
                cur_mouse_vals = vessel_long_vector(rows,7); %store for long format
        end
        
        %plot vector and patch
        ul = mean_vec+sem_vec;
        ll = mean_vec-sem_vec;
        patch([time_vec fliplr(time_vec)],[ul fliplr(ll)],'k','FaceColor',colors(c,:),'EdgeColor','none','FaceAlpha',alphas(c));
        plot(time_vec,mean_vec,[linetypes{c} 'k'],'Color',colors(c,:),'LineWidth',1);

    end
    %add figure elements
    xlim([-3 15])
    ylim(vector_ylimits(dt,:))
    xticks(0:5:15)
    box off
    set(gca,'FontSize',font_size)
    set(gca,'FontName','Arial')
    xlabel('time (s)','FontSize',font_size)
    ylabel(ylabels{dt},'FontSize',font_size)

    %plot response analysis ranges
    ylineloc = vector_ylimits(dt,1) + 0.05*(vector_ylimits(dt,2)-vector_ylimits(dt,1));
    ylineh = 0.02*(vector_ylimits(dt,2)-vector_ylimits(dt,1));
    switch dt
        case 1
            plot([0.5 6],[ylineloc ylineloc],'k')
            plot([0.5 0.5],ylineloc+[ylineh -ylineh],'k')
            plot([6 6],ylineloc+[ylineh -ylineh],'k')
            text(3.25,ylineloc,'response','FontSize',font_size,'FontName','Arial','HorizontalAlignment','center','VerticalAlignment','bottom')
        case 2
            plot([0.5 6],[ylineloc ylineloc],'k')
            plot([0.5 0.5],ylineloc+[ylineh -ylineh],'k')
            plot([6 6],ylineloc+[ylineh -ylineh],'k')
            text(3.25,ylineloc,'response','FontSize',font_size,'FontName','Arial','HorizontalAlignment','center','VerticalAlignment','bottom')
        case 3
            plot([0.5 2.5],[ylineloc ylineloc],'k')
            plot([2.5 12],[ylineloc ylineloc],'k')
            plot([0.5 0.5],ylineloc+[ylineh -ylineh],'k')
            plot([2.5 2.5],ylineloc+[ylineh -ylineh],'k')
            plot([12 12],ylineloc+[ylineh -ylineh],'k')
            text(1.5,ylineloc,'rise','FontSize',font_size,'FontName','Arial','HorizontalAlignment','center','VerticalAlignment','bottom')
            text(7.25,ylineloc,'recovery','FontSize',font_size,'FontName','Arial','HorizontalAlignment','center','VerticalAlignment','bottom')
        case 4
            plot([0.5 6],[ylineloc ylineloc],'k')
            plot([0.5 0.5],ylineloc+[ylineh -ylineh],'k')
            plot([6 6],ylineloc+[ylineh -ylineh],'k')
            text(3.25,ylineloc,'response','FontSize',font_size,'FontName','Arial','HorizontalAlignment','center','VerticalAlignment','bottom')
    end
    print(gcf,'-vector','-dsvg',['C:\Users\misaa\Desktop\2p vector ' dtypenames{dt} '.svg'])
end

% save(fullfile(base_folder,'vel_avgs_2p.mat'),'vel_means')


%% (Plotting) create delta delta vector and delta delta box plots

%set time ranges for analysis of neuron/vessel response
neuron_response_inds = find(time_vec>0.5 & time_vec<=6);
vessel_response_inds = find(time_vec>0.5 & time_vec<=12);
vessel_rise_inds = find(time_vec>0.5 & time_vec<=2.5);
vessel_recovery_inds = find(time_vec>2.5 & time_vec<=12);

do_stats = false;
do_plots = true;

vector_ylimits = [-7 7; -7 7; -10 10; -1 1.5]; %[pos neuron; neg neuron; velocity; width]
% box_ylimits = [-30 50; -30 50; -25 25; -5 7];
box_ylimits = [-30 50; -30 50; -18 13; -5 7];
lme_delta_ps = nan(5,1);
cur_comparison = 0;
dtypenames = {'posneurons','negneurons', 'velocities','widths'};
ylabels_deltas = {'mean % \DeltaF/F', 'mean % \DeltaF/F', 'mean rise % \DeltaV/V', 'mean recovery % \DeltaV/V', 'mean % \DeltaW/W'};
curlabel = 0;
conds_to_compare = [2 4];  %control pre/post, psilo pre/post
% conds_to_compare = [2 8];  %control pre/post, psilo+mdl pre/post
% conds_to_compare = [4 8];  %psilo pre/post, psilo+mdl pre/post
cur_dt = 0;
for dt = [1 2 3 4] %posneuron, negneuron, velocity, width
    %allocate space for delta values
    pre_vals = nan(1,8,4);
    post_vals = nan(1,8,4);
    delta_vals = nan(1,8,4);
    Data = [];
    Treatment = [];
    Mouse = [];
    Timepoint = [];
    Cell = [];
    if do_plots
        %create figure for vectors
        figure('Position',[100 100 vector_plot_size])
        plot([0 20],[0 0],'--k')
        hold on
        patch([0 6 6 0],[-100 -100 100 100],'k','EdgeColor','none','FaceAlpha',0.1)
    end

    for c = conds_to_compare
        switch dt
            case 1
                rows = find(pre_pos_neurons & neuron_long_vector(:,8,1)==c); %get all "either_pos" neurons for this condition
                nr = length(rows);
                cur_cells = neuron_long_vector(rows,8,1);
                cur_mouse_vals = neuron_long_vector(rows,5,1); %store for long format
                deltas = neuron_long_vector(rows,2,:) - neuron_long_vector(rows-1,2,:); %delta delta
                mean_vec = squeeze(mean(deltas,1,'omitnan'))'; %mean dF/F vector
                std_vec = squeeze(std(deltas,0,1,'omitnan'))';
                sem_vec = std_vec/sqrt(nr);
                cur_post_vals = squeeze(mean(neuron_long_vector(rows,2,neuron_response_inds),3,'omitnan')); %delta mean dF/F during stim
                cur_pre_vals = squeeze(mean(neuron_long_vector(rows-1,2,neuron_response_inds),3,'omitnan'));
                cur_delta_vals = cur_post_vals - cur_pre_vals;  %delta min dV/V during/after stim

            case 2
                rows = find(pre_neg_neurons & neuron_long_vector(:,8,1)==c); %get all "either_pos" neurons for this condition
                nr = length(rows);
                cur_cells = neuron_long_vector(rows,8,1);
                cur_mouse_vals = neuron_long_vector(rows,5,1); %store for long format
                deltas = neuron_long_vector(rows,2,:) - neuron_long_vector(rows-1,2,:); %delta delta
                mean_vec = squeeze(mean(deltas,1,'omitnan'))'; %mean dF/F vector
                std_vec = squeeze(std(deltas,0,1,'omitnan'))';
                sem_vec = std_vec/sqrt(nr);
                cur_post_vals = squeeze(mean(neuron_long_vector(rows,2,neuron_response_inds),3,'omitnan')); %delta mean dF/F during stim
                cur_pre_vals = squeeze(mean(neuron_long_vector(rows-1,2,neuron_response_inds),3,'omitnan'));
                cur_delta_vals = cur_post_vals - cur_pre_vals;  %delta min dV/V during/after stim

            case 3
                rows = find(all_vels & vessel_long_vector(:,10,1)==c); %get all velocity-measured vessels for this condition
                nr = length(rows);
                cur_cells = vessel_long_vector(rows,8,1);
                cur_mouse_vals = vessel_long_vector(rows,7,1); %store for long format
                deltas = vessel_long_vector(rows,2,:) - vessel_long_vector(rows-1,2,:);
                mean_vec = squeeze(mean(deltas,1,'omitnan'))'; %mean dF/F vector
                std_vec = squeeze(std(deltas,0,1,'omitnan'))';
                sem_vec = std_vec/sqrt(nr);
                %get mean of peak/undershoot regions by finding mean during stim region and mean during post stim region
                cur_post_vals = squeeze(mean(vessel_long_vector(rows,2,vessel_rise_inds),3,'omitnan'));
                cur_pre_vals = squeeze(mean(vessel_long_vector(rows-1,2,vessel_rise_inds),3,'omitnan'));
                cur_delta_vals = cur_post_vals - cur_pre_vals;  %delta max dV/V during/after stim
                cur_post_vals(:,:,2) = squeeze(mean(vessel_long_vector(rows,2,vessel_recovery_inds),3,'omitnan'));
                cur_pre_vals(:,:,2) = squeeze(mean(vessel_long_vector(rows-1,2,vessel_recovery_inds),3,'omitnan'));
                cur_delta_vals(:,:,2) = cur_post_vals(:,:,2) - cur_pre_vals(:,:,2);  %delta min dV/V during/after stim

            case 4
                rows = find(all_widths & vessel_long_vector(:,10,1)==c); %get all width-measured vessels for this condition
                nr = length(rows);
                cur_cells = vessel_long_vector(rows,8,1);
                cur_mouse_vals = vessel_long_vector(rows,7,1); %store for long format
                deltas = vessel_long_vector(rows,5,:) - vessel_long_vector(rows-1,5,:);
                mean_vec = squeeze(mean(deltas,1,'omitnan'))'; %mean dF/F vector
                std_vec = squeeze(std(deltas,0,1,'omitnan'))';
                sem_vec = std_vec/sqrt(nr);
                cur_post_vals = squeeze(mean(vessel_long_vector(rows,5,vessel_response_inds),3,'omitnan'));
                cur_pre_vals = squeeze(mean(vessel_long_vector(rows-1,5,vessel_response_inds),3,'omitnan'));
                cur_delta_vals = cur_post_vals - cur_pre_vals;  %delta min dV/V during/after stim

        end
        if do_plots
            %plot vector and patch
            ul = mean_vec+sem_vec;
            ll = mean_vec-sem_vec;
            patch([time_vec fliplr(time_vec)],[ul fliplr(ll)],'k','FaceColor',colors(c,:),'EdgeColor','none','FaceAlpha',alphas(c));
            plot(time_vec,mean_vec,[linetypes{c} 'k'],'Color',colors(c,:),'LineWidth',1);
        end

        %store deltas
        if nr>size(delta_vals,1)
            pre_vals(end+1:nr,:,:) = nan;
            post_vals(end+1:nr,:,:) = nan;
            delta_vals(end+1:nr,:,:) = nan;
        end
        pre_vals(1:nr,c,1) = cur_pre_vals(:,:,1);
        post_vals(1:nr,c,1) = cur_post_vals(:,:,1);
        delta_vals(1:nr,c,1) = cur_delta_vals(:,:,1);
        if dt==3 %flow recovery
            pre_vals(1:nr,c,2) = cur_pre_vals(:,:,2);
            post_vals(1:nr,c,2) = cur_post_vals(:,:,2);
            delta_vals(1:nr,c,2) = cur_delta_vals(:,:,2);
        end

        %store measurement and metadata in long format for stats
        Data = [Data; cur_pre_vals];
        Mouse = [Mouse; cur_mouse_vals];
        Treatment = [Treatment; c*ones(nr,1)];
        Timepoint = [Timepoint; ones(nr,1)];
        Cell = [Cell; cur_cells];

        Data = [Data; cur_post_vals];
        Mouse = [Mouse; cur_mouse_vals];
        Treatment = [Treatment; c*ones(nr,1)];
        Timepoint = [Timepoint; 2*ones(nr,1)];
        Cell = [Cell; cur_cells];
    end
    if do_plots
        %add figure elements
        xlim([-3 15])
        ylim(vector_ylimits(dt,:))
        xticks(0:5:15)
        box off
        set(gca,'FontSize',font_size)
        set(gca,'FontName','Arial')
        xlabel('time (s)','FontSize',font_size)
        ylabel(ylabels_deltas{dt},'FontSize',font_size)

        % print(gcf,'-vector','-dsvg',fullfile(base_folder,['2p deltavector' dtypenames{dt} ' ' num2str(i) '.svg']))
    end

    %create figure for pre/post/delta boxplots
    if dt==3
        nn = 2;
    else 
        nn = 1;
    end
    boxplot_vals = [pre_vals(:,conds_to_compare(1),:) post_vals(:,conds_to_compare(1),:) delta_vals(:,conds_to_compare(1),:) nan(size(delta_vals(:,conds_to_compare(1),:))) pre_vals(:,conds_to_compare(2),:) post_vals(:,conds_to_compare(2),:) delta_vals(:,conds_to_compare(2),:)];
    boxplot_colors = [colors(conds_to_compare(1)-1,:); colors(conds_to_compare(1),:); colors(conds_to_compare(1),:); colors(conds_to_compare(2)-1,:); colors(conds_to_compare(2),:); colors(conds_to_compare(2),:)];
    for i = 1:nn
        curlabel = curlabel+1;

        if do_plots
            figure('Position',[100 100 box_plot_size])
            plot([0 8],[0 0],':k')
            hold on
            x = kron([1 2 3 5 6 7],ones(size(boxplot_vals,1),1));
            for j = 1:6
                swarmchart(x(:,j),boxplot_vals(:,x(1,j),i),1,0.8*boxplot_colors(j,:),'filled','o','XJitterWidth',0.25,'MarkerFaceAlpha',0.7,'MarkerEdgeColor','none')
            end
            bp = boxplot(boxplot_vals(:,:,i),'Symbol','');
            set(bp(6,:), 'Color', 'k'); %median
            set(bp(5,:), 'Color', 'k'); %box
            patch([3.75 4.25 4.25 3.75],[box_ylimits(dt,1)*1.1 box_ylimits(dt,:)*1.1 box_ylimits(dt,2)*1.1],'w','EdgeColor','None')
            plot([3 7],1.05*[box_ylimits(dt,2) box_ylimits(dt,2)],'k')
            plot([3 3],[0.95*box_ylimits(dt,2) 1.05*box_ylimits(dt,2)],'k')
            plot([7 7],[0.95*box_ylimits(dt,2) 1.05*box_ylimits(dt,2)],'k')
            text(3,1.1*box_ylimits(dt,1),'\Delta','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','center')
            text(7,1.1*box_ylimits(dt,1),'\Delta','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','center')
            ylim(box_ylimits(dt,:))
            xlim([0 8])
            set(gca,'FontSize',font_size)
            set(gca,'FontName','Arial')
            ylabel(ylabels_deltas{curlabel},'FontSize',font_size)
            % xticks(1:3)
            % yticks(-50:10:50)
            xticklabels({'pre','post','','','pre','post',''})
            % xticklabels({'p','p','d','','p','p','d'})
            ax = gca;
            ax.Clipping = "off";
            box off
        end

        if do_stats
            %run LME for this dataset
            cur_comparison = cur_comparison+1;
            lmedata.measurement = Data(:,:,i);
            lmedata.mouse = categorical(Mouse);
            lmedata.treatment = categorical(Treatment);
            lmedata.timepoint = categorical(Timepoint);
            lmedata.cell = categorical(Cell);
            lmetable = struct2table(lmedata);
            lme_model = fitlme(lmetable, 'measurement ~ treatment*timepoint + (1|mouse) + (1|cell:mouse)');

            % ANOVA p-value for treatment effect
            % anova_results = anova(lme_model); 
            lme_delta_ps(cur_comparison,:) = lme_model.Coefficients.pValue(4);


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

        if do_plots
            text(5,1.1*box_ylimits(dt,2),'n.s. (p=0.00)','FontSize',font_size,'VerticalAlignment','bottom','HorizontalAlignment','center')

            % print(gcf,'-vector','-dsvg',fullfile(base_folder,['2p boxplot ' dtypenames{dt} ' ' num2str(i) '.svg']))
        end
        

    end
end
if do_stats
    all_ps = [lme_baseline_ps; lme_delta_ps]
    [~, ~, ~, adj_b_ps] = fdr_bh(all_ps)
end


%% (Plotting) create bar chart of classifications
barnames = {'neurons', 'velocities','widths'};
bar_ylims = [100 100 100];
% conds_to_compare = [1 2 3 4]; %control pre/post, psilo pre/post
conds_to_compare = [7 8];
for dt = [1 2 3] %neuron, velocity, width
    figure('Position',[100 100 box_plot_size])
    classes = nan(8,3); %[cond, class]
    for c = conds_to_compare 
        switch dt
            case 1
                classes(c,1) = sum(pos_neurons & neuron_long_vector(:,8,1)==c); %get all "POS" neurons for this condition
                classes(c,2) = sum(neg_neurons & neuron_long_vector(:,8,1)==c); %get all "NEG" neurons for this condition
                classes(c,3) = sum(none_neurons & neuron_long_vector(:,8,1)==c); %get all "NONE" neurons for this condition
            case 2 
                classes(c,1) = sum(pos_vels & vessel_long_vector(:,10,1)==c); %get all "POS" vessels for this condition
                classes(c,2) = sum(neg_vels & vessel_long_vector(:,10,1)==c); %get all "NEG" vessels for this condition
                classes(c,3) = sum(none_vels & vessel_long_vector(:,10,1)==c); %get all "NONE" vessels for this condition
            case 3
                classes(c,1) = sum(pos_widths & vessel_long_vector(:,10,1)==c); %get all "POS" vessels for this condition
                classes(c,2) = sum(neg_widths & vessel_long_vector(:,10,1)==c); %get all "NEG" vessels for this condition
                classes(c,3) = sum(none_widths & vessel_long_vector(:,10,1)==c); %get all "NONE" vessels for this condition
        end
        
        num_classes = sum(~isnan(classes(c,:)));
        for i = 1:num_classes
            [phat, pci] = binofit(classes(c,i), sum(classes(c,1:num_classes)));  
            phats(c,i) = 100*phat;
            ci_lower(c,i) = 100*(pci(1)-phat); % lower error
            ci_upper(c,i) = 100*(phat-pci(2)); % upper error
        end

    end
    classes_prct = 100*classes./repmat(sum(classes,2,'omitnan'),[1 3]);
    classes_prct = classes_prct(conds_to_compare,:);
    b = bar(classes_prct');
    hold on
    ylim([0 bar_ylims(dt)])
    phats = phats(conds_to_compare,:);
    ci_lower = ci_lower(conds_to_compare,:);
    ci_upper = ci_upper(conds_to_compare,:);
    % x = [1; 2; 3; 4] + [-0.27 -0.09 0.09 0.27];
    x = [1; 2; 3; 4] + [-0.14 0.14];
    e = errorbar(x', phats, ci_lower, ci_upper, 'k.', 'LineWidth', 0.5);
    for i = 1:length(e)
        e(i).CapSize = 1;
    end
    for co = 1:length(conds_to_compare)
        b(co).FaceColor = colors(conds_to_compare(co),:);
        b(co).FaceAlpha = 0.7;
    end
    set(gca,'FontSize',font_size)
    set(gca,'FontName','Arial')
    xticklabels({'pos.','neg.','none'})
    xlim([0.5 3.5])
    if dt==1
        ylabel('% of neuron responsivities ','FontSize',font_size)
    else
        ylabel('% of vessel responsivities','FontSize',font_size)
    end
    box off

    % print(gcf,'-vector','-dsvg',fullfile(base_folder,['2p classes ' barnames{dt} '.svg']))
end




%% load data to create methods/example panels
%find example raw for WT condition that has many responding cells and vessels
% num_rn = sum(responding_neurons_on,3);
% num_rv = sum(responding_vels,3);
% best_idx = (num_rv==3 & num_rn==6); %examples with all 3 vessels and all 6 neurons positively responding
% best_wt_mouse = find(best_idx(1,:),1);
% best_dir = fileparts(velnames{1,best_wt_mouse,1});
% best_raw_data_name = fullfile(best_dir,ls(fullfile(best_dir,'*.pmt.dat')));
% best_trig_name = fullfile(best_dir,'Ch1_triggers_extracted.mat');

%find example data that Rick has linescan path data for
best_dir = fileparts(velnames{4,4,1});
best_raw_data_name = fullfile(best_dir,ls(fullfile(best_dir,'*.pmt.dat')));
best_trig_name = fullfile(best_dir,'Ch1_triggers_extracted.mat');


% metadata
numChannels = 4;
samplesPerFrame = 15000;
linestoinclude = 4000;
channel = 4;

%load the data (not the whole thing, since the file is too large)
datapoints = samplesPerFrame*linestoinclude*numChannels;
fileID = fopen(best_raw_data_name);
data = fread(fileID,datapoints,'int16');

%format data for images
linescanData = reshape(data,numChannels,samplesPerFrame,[]); %reshape for: %[channels, samples, lines]

vesselData = permute(linescanData(4,:,1:4000),[2 3 1]); %reshape for [samples, lines, channels]
vesselData = vesselData/3500;
vesselData(vesselData<0) = 0;
vesselData(vesselData>1) = 1;
vesselData = imresize(vesselData,[1080 1920]);

neuralData = permute(linescanData(2,:,1:4000),[2 3 1]); 
neuralData = neuralData/4000;
neuralData(neuralData<0) = 0;
neuralData(neuralData>1) = 1;
neuralData = imresize(neuralData,[1080 1920]);

stimData = permute(linescanData(1,:,1:4000),[2 3 1]); 
stimData = stimData/-30000;
stimData(stimData<0) = 0;
stimData(stimData>1) = 1;
stimData = imresize(stimData,[1080 1920]);

%combine vessel (red) and neural (green) channels
combImage = repmat(vesselData*2,[1 1 3]);
combImage(:,:,1) = combImage(:,:,1)*colors(8,1);
combImage(:,:,2) = combImage(:,:,2)*colors(8,2);
combImage(:,:,3) = combImage(:,:,3)*colors(8,3);
combImage(:,:,2) = combImage(:,:,2) + neuralData*2;


%load example traces
%flow:
curfilename = velnames{4,4,3};
load(curfilename)
exp_folder = fileparts(curfilename);
curtrigfilename = ls(fullfile(exp_folder,'Ch1_triggers_extracted*'));
assert(size(curtrigfilename,1)==1,'unexpected # of trigger files')
curtrigfilename = fullfile(exp_folder,curtrigfilename);
load(curtrigfilename)
channelAvg = medfilt1(channelAvg,500);
minTrig = min(channelAvg);
maxTrig = max(channelAvg);
channelAvg = channelAvg>mean([minTrig maxTrig]);
risefall = find(diff(channelAvg)~=0);
trig_inds = round(risefall(1:2:end)/50);
num_reps = length(risefall)/2;
assert((num_reps==25)|(num_reps==20),'unexpected number of stimuli')
velocities = Result(:,3);
medvel = median(velocities);
if medvel<0
    medvel = -medvel;
    velocities = -velocities;
end
velocities(velocities>(1.6*medvel)) = nan;
velocities(velocities<(0.5*medvel)) = nan;
for i = 1:3
    velocities_smooth = movmedian(velocities,20);
    bad_idx = abs(velocities_smooth-velocities)>(0.2*velocities_smooth);
    velocities(bad_idx) = nan;
    velocities = fillmissing(velocities,"linear",'MissingLocations',isnan(velocities));
end
velocities_smooth = movmean(velocities,3);
tmps = nan(5,length(time_vec));
for r = 1:5
    inds = (trig_inds(r)-10):(trig_inds(r)+num_pts-11);
    inds(inds>length(velocities_smooth)) = [];
    tmp = velocities_smooth(inds);
    tmps(r,1:length(tmp)) = tmp;
end
tmps = tmps./repmat(mean(tmps(:,baseline_inds),2,'omitnan'),[1 num_pts]);
tmps = 100*(tmps-1);

%diameter:
curfilename = diamnames{4,4,3};
load(curfilename)
exp_folder = fileparts(curfilename);
curtrigfilename = ls(fullfile(exp_folder,'Ch1_triggers_extracted*'));
assert(size(curtrigfilename,1)==1,'unexpected # of trigger files')
curtrigfilename = fullfile(exp_folder,curtrigfilename);
load(curtrigfilename)
channelAvg = medfilt1(channelAvg,500);
minTrig = min(channelAvg);
maxTrig = max(channelAvg);
channelAvg = channelAvg>mean([minTrig maxTrig]);
risefall = find(diff(channelAvg)~=0);
trig_inds = round(risefall(1:2:end)/50);
num_reps = length(risefall)/2;
assert((num_reps==25)|(num_reps==20),'unexpected number of stimuli')
widths = results.vesselDiameterInMicronsLineVec;
medwidth = median(widths);
if medwidth<0
    medwidth = -medwidth;
    widths = -widths;
end
widths(widths>(1.6*medwidth)) = nan;
widths(widths<(0.5*medwidth)) = nan;
for i = 1:3
    widths_smooth = movmedian(widths,20);
    bad_idx = abs(widths_smooth-widths)>(0.2*widths_smooth);
    widths(bad_idx) = nan;
    widths = fillmissing(widths,"linear",'MissingLocations',isnan(widths));
end
widths_smooth = movmean(widths,30);
d_tmps = nan(5,length(time_vec));
for r = 1:5
    inds = (trig_inds(r)-10):(trig_inds(r)+num_pts-11);
    inds(inds>length(widths_smooth)) = [];
    tmp = widths_smooth(inds);
    d_tmps(r,1:length(inds)) = tmp;
end
d_tmps = d_tmps./repmat(mean(d_tmps(:,baseline_inds),2,'omitnan'),[1 num_pts]);
d_tmps = 100*(d_tmps-1);

%neurons:
curfilename = neuronnames{4,4,4};
cur_image = read_file(curfilename);
[h,w,nf] = size(cur_image);
cur_image(cur_image==0) = nan;
fluo_vec_1 = mean(cur_image,2,'omitnan');
fluo_vec_1 = reshape(fluo_vec_1,[1, numel(fluo_vec_1)]);
fluo_vec_1 = imresize(fluo_vec_1,[1 length(velocities_smooth)],"bicubic");
curtrigfilename = ls(fullfile(exp_folder,'Ch1_triggers_extracted*'));
if ~isempty(curtrigfilename)
    assert(size(curtrigfilename,1)==1,'unexpected # of trigger files')
    curtrigfilename = fullfile(exp_folder,curtrigfilename);
    load(curtrigfilename)
    channelAvg = medfilt1(channelAvg,500);
    minTrig = min(channelAvg);
    maxTrig = max(channelAvg);
    channelAvg = channelAvg>mean([minTrig maxTrig]);
    risefall = find(diff(channelAvg)~=0);
    trig_inds = round(risefall(1:2:end)/50);
    num_reps = length(risefall)/2;
    assert((num_reps==25)|(num_reps==20),'unexpected number of stimuli')
end
n_tmps = nan(5,length(time_vec));
for r = 1:5
    inds = (trig_inds(r)-10):(trig_inds(r)+num_pts-11);
    inds(inds>length(fluo_vec_1)) = [];
    tmp = fluo_vec_1(inds);
    n_tmps(r,1:length(inds)) = tmp;
end
n_tmps = n_tmps./repmat(mean(n_tmps(:,baseline_inds),2,'omitnan'),[1 num_pts]);
n_tmps = 100*(n_tmps-1);


%% (Plotting) plot example line scan data
%(panels A) experimental setup
axes('Position',[0.04 0.81 0.17 0.17])
text(0,0,'experimental setup','HorizontalAlignment','center','FontSize',6)
xlim([-1 1])
ylim([-1 1])
box off
axis off
set(gca,'FontSize',6)

%(panels B) linescan path
axes('Position',[0.24 0.81 0.17 0.17])
text(0,0,'example linescan path','HorizontalAlignment','center','FontSize',6)
xlim([-1 1])
ylim([-1 1])
box off
axis off
set(gca,'FontSize',6)


%(panel C) example linescan raw data
stimVec = medfilt1(mean(stimData),100);
axes('Position',[0.5 0.81 0.26 0.17])
imshow(combImage);
hold on
plot(1:1920, -75*stimVec-50,'k')
axis on
% xticks(1:167:167*21.1)
% xticklabels(0:20)
xlabel('time (s)')
ylabel('linescan (pixels)')
% text(167*12,25,'vessel velocity 1','Color','k','FontSize',6)
% text(167*12,115,'neuron 1','Color','k','FontSize',6)
% text(167*12,205,'neuron 2','Color','k','FontSize',6)
% text(167*12,295,'vessel width 1','Color','k','FontSize',6)
% text(167*12,395,'vessel velocity 2','Color','k','FontSize',6)
% text(167*12,485,'neuron 3','Color','k','FontSize',6)
% text(167*12,575,'neuron 4','Color','k','FontSize',6)
% text(167*12,665,'vessel width 2','Color','k','FontSize',6)
% text(167*12,765,'vessel velocity 3','Color','k','FontSize',6)
% text(167*12,855,'neuron 5','Color','k','FontSize',6)
% text(167*12,945,'neuron 6','Color','k','FontSize',6)
% text(167*12,1035,'vessel width 3','Color','k','FontSize',6)
set(gca,'FontSize',6)
box off
%text(0,0,'example linescan raw data','HorizontalAlignment','center','FontSize',6)
xlim([1 1920])
ylim([-125 1080])
box off
axis off
set(gca,'FontSize',6)


%(panel D) example extracted raw data 
axes('Position',[0.80 0.81 0.15 0.17])
x_new = linspace(min(time_vec), max(time_vec), 300);
stimVec2 = double(x_new>=0 & x_new<=6);
plot(x_new, 25*stimVec2+350,'k')
hold on
plot(time_vec,200*ones(size(time_vec)),'--k','Color',colors(1,:))
plot(time_vec,100*ones(size(time_vec)),'--k','Color',colors(1,:))
plot(time_vec,zeros(size(time_vec)),'--k','Color',colors(1,:))
plot(time_vec,200+n_tmps'/3,'Color',[0 0.8 0])
plot(time_vec,100+tmps','Color',colors(8,:))
plot(time_vec,000+d_tmps','Color',colors(8,:))
xlim([-3 15])
ylim([-25 375])
xlabel('time (s)')
box off
ax = gca;
set(ax, 'YColor', 'none');
set(ax, 'YTick', []);
set(ax, 'YLabel', []);
set(ax, 'Clipping', 'off');
plot([-4 -4],[0 50],'k')
plot([-4 -4],[100 150],'k')
plot([-4 -4],[200 250],'k')
% axis off
%text(0,300,'example traces','HorizontalAlignment','center','FontSize',6)
set(gca,'FontSize',6)

% print(gcf,'-vector','-dsvg',fullfile(base_folder,'2p examples.svg'))
    


%% unused extras

% deconvolved impulse response function
stim = zeros(num_pts,1);
stim(time_vec>=0 & time_vec<=6) = 1;
figure('Position',[100 100 1400 250])


%plot input
subplot(1,4,1)
plot(time_vec,zeros(size(time_vec)),'--k')
hold on
plot(time_vec,stim,'k','LineWidth',1)
ylim([-0.2 1.2])
title('input (stimulus)')
stim(stim==0) = [];
stim = stim/sum(stim);


%plot output
subplot(1,4,2)
exclude_idx = repmat((inclusive_positive_vessels~=1),[1 1 1 num_pts]);
dvv_avg = dvv; %[num_conds, num_mice, num_vessels, numpts]
dvv_avg(exclude_idx) = nan;
dvv_avg = squeeze(mean(dvv_avg,3,'omitnan')); %[num_conds, num_mice, numpts]
num_plotted_mice = sum(~isnan(dvv_avg(:,:,1)),2)';
avg1 = squeeze(mean(dvv_avg,2,'omitnan')); %[num_conds, numpts]
std1 = squeeze(std(dvv_avg,[],2,'omitnan'))';
sem1 = std1./sqrt(repmat(num_plotted_mice,[num_pts 1]));
plot([0 20],[0 0],'--k')
hold on
patch([0 6 6 0],[-100 -100 100 100],'k','EdgeColor','none','FaceAlpha',0.15)
for c = 3:4
    x = time_vec;
    y = avg1(c,:);
    ul = y + sem1(:,c)';
    ll = y - sem1(:,c)';
    patch([x fliplr(x)],[ul fliplr(ll)],[linetypes{c} 'k'],'FaceColor',colors(c,:),'EdgeColor','none','FaceAlpha',alphas(c));
    plot(x,y,[linetypes{c} 'k'],'Color',colors(c,:),'LineWidth',1);
end
xlim([-3 15])
ylim([-9 13])
xticks(0:5:15)
xlabel('time (s)')
ylabel('\DeltaV/V')
title('output (flow)')
legend({'','','','Pre-Psilocybin','','Post-Psilocybin'})


%plot deconvolution
deconv_colors = colors;
deconv_colors(3,:) = [1 0.5 1];
deconv_colors(4,:) = [0.9 0 0.9];
subplot(1,4,3)
hrf_max_time = 6;
hrf_max_pts = round(hrf_max_time/P);
hrfs = nan(num_conds,num_mice,hrf_max_pts);
for c = 3:4
    for m = 1:num_mice
        y = squeeze(dvv_avg(c,m,:));
        if ~all(isnan(y))
            [x,r] = deconv(y,stim);
            x(1:10) = [];
            x = x(1:hrf_max_pts);
            hrfs(c,m,:) = movmean(x,3);
        end
    end
end
tv = 0:P:(length(x)-1)*P;
hrfs_avg = squeeze(mean(hrfs,2,'omitnan'));
hrfs_std = squeeze(std(hrfs,[],2,'omitnan'));
hrfs_sem = hrfs_std./sqrt(repmat(num_plotted_mice',[1 length(x)]));
plot(0:6,zeros(1,7),'--k')
hold on
plot([0 0],[-100 100],'k','LineWidth',1)
for c = 3:4
    ul = hrfs_avg(c,:) + hrfs_sem(c,:);
    ll = hrfs_avg(c,:) - hrfs_sem(c,:);
    patch([tv fliplr(tv)],[ll fliplr(ul)],[linetypes{c} 'k'],'EdgeColor','none','FaceColor',deconv_colors(c,:),'FaceAlpha',alphas(c))
    plot(tv,hrfs_avg(c,:),[linetypes{c} 'k'],'Color',deconv_colors(c,:),'LineWidth',1)
end
plot([0 0],[0 1],'k','LineWidth',2)
ylim([-25 55])
xlim([-1 hrf_max_time+1])
xlabel('time (s)')
title('HRF = deconv(output,input)')
legend({'','','','pre-psilocybin HRF','','post-psilocybin HRF'})


%convolve to get original signals back
subplot(1,4,4)
plot([0 20],[0 0],'--k')
hold on
patch([0 6 6 0],[-100 -100 100 100],'k','EdgeColor','none','FaceAlpha',0.15)
for c = 3:4
    x = time_vec;
    y = avg1(c,:);
    ul = y + sem1(:,c)';
    ll = y - sem1(:,c)';
    patch([x fliplr(x)],[ul fliplr(ll)],[linetypes{c} 'k'],'FaceColor',colors(c,:),'EdgeColor','none','FaceAlpha',alphas(c));
    plot(x,y,[linetypes{c} 'k'],'Color',colors(c,:),'LineWidth',1);
end
for c = 3:4
    y = squeeze(mean(hrfs(c,:,:),2,'omitnan'));
    y = conv(y,stim);
    x = time_vec;
    x(x<0) = [];
    x = x(1:length(y));
    plot(x,y,[linetypes{c} 'k'],'Color',deconv_colors(c,:),'LineWidth',1);
end
xlim([-3 15])
ylim([-9 13])
xticks(0:5:15)
xlabel('time (s)')
ylabel('\DeltaV/V')
title('simulated = conv(hrf,input)')
legend({'','','','Pre-Psilocybin','','Post-Psilocybin','Pre-Psilocybin (simulated)','Post-Psilocybin (simulated)'})



%(panel) plot vessel flow deltas for all positively responding vessels
plot(time_vec,zeros(size(time_vec)),'--k')
hold on
patch([0 6 6 0],[-100 -100 100 100],'k','EdgeColor','none','FaceAlpha',0.1)
num_plotted_mice = nan(1,4);
for c = 2:2:4
    d_avg = dvv;

    %exclude all non-positive responders
    nonincluded_idx = repmat(inclusive_positive_vessels,[1 1 1 num_pts])~=1;
    d_avg(nonincluded_idx) = nan;
    d_avg_delta = d_avg(c,:,:,:) - d_avg(c-1,:,:,:);  %[1 mice vessel trace]

    %averaging all deltas within each mouse
    d_avg_delta = permute(d_avg_delta,[4 1 2 3]); %[trace, 1, mouse, vessel]
    d_avg_delta = squeeze(mean(d_avg_delta,4,'omitnan')); %[trace, mouse]

    num_plotted_mice(c) = sum(~isnan(d_avg_delta(1,:)));
    d_avg_mean = mean(d_avg_delta,2,'omitnan')';
    d_avg_std = std(d_avg_delta,0,2,'omitnan')';
    d_avg_sem = d_avg_std/sqrt(num_plotted_mice(c));

    ul = d_avg_mean + d_avg_sem;
    ll = d_avg_mean - d_avg_sem;
    patch([time_vec fliplr(time_vec)],[ul fliplr(ll)],'k','FaceColor',colors(c,:),'EdgeColor','none','FaceAlpha',0.15)
    plot(time_vec,d_avg_mean,'Color',colors(c,:),'LineWidth',1)
end
xlim([-3 15])
ylim([-5 10])
xlabel('time (s)')
ylabel('\DeltaV/V (%)')
box off
set(gca,'FontSize',6)
title('velocity \Delta')


%neuron deltas show too much variability in pre/post timepoints, even control
%treatment shows changes
%(panel) neuron deltas
plot(time_vec,zeros(size(time_vec)),'--k')
hold on
patch([0 6 6 0],[-100 -100 100 100],'k','EdgeColor','none','FaceAlpha',0.1)
num_plotted_mice = nan(1,4);
for c = 2:2:4
    d_avg = dff;

    %exclude all non-positive responders
    nonincluded_idx = repmat(inclusive_positive_neurons,[1 1 1 num_pts])~=1;
    d_avg(nonincluded_idx) = nan;
    d_avg_delta = d_avg(c,:,:,:) - d_avg(c-1,:,:,:);  %[1 mice neuron trace]

    %averaging all deltas, within each mouse
    d_avg_delta = permute(d_avg_delta,[4 1 2 3]); %[trace, 1, mouse, neuron]
    d_avg_delta = squeeze(mean(d_avg_delta,4,'omitnan')); %[trace, mouse]

    num_plotted_mice(c) = sum(~isnan(d_avg_delta(1,:)));
    d_avg_mean = mean(d_avg_delta,2,'omitnan')';
    d_avg_std = std(d_avg_delta,0,2,'omitnan')';
    d_avg_sem = d_avg_std/sqrt(num_plotted_mice(c));

    ul = d_avg_mean + d_avg_sem;
    ll = d_avg_mean - d_avg_sem;
    patch([time_vec fliplr(time_vec)],[ul fliplr(ll)],'k','FaceColor',colors(c,:),'EdgeColor','none','FaceAlpha',0.15)
    plot(time_vec,d_avg_mean,'Color',colors(c,:),'LineWidth',1)
end
xlim([-3 15])
ylim([-15 15])
xlabel('time (s)')
ylabel('\DeltaF/F (%)')
box off
set(gca,'FontSize',6)
title('neural activity \Delta')


