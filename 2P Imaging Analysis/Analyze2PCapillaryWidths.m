%% set base directory
base_folder = 'D:\2P Data\New file organization\';

%%

% --- Parameters
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
SN_colors = SN_colors / 255;

colors = [SN_colors(3,:);...% pre control - grey
          SN_colors(1,:);...% post control - black
          SN_colors(9,:);...% pre psilo - pale-blue
          SN_colors(10,:);... % post psilo - blue
          SN_colors(6,:);...% pre mdl+psilo - pale-orange
          SN_colors(5,:)];   % post mdl+psilo - red

mouseDirs = dir(base_folder);
mouseDirs = mouseDirs([mouseDirs.isdir] & ~startsWith({mouseDirs.name}, '.') & ~cellfun(@isempty, regexp({mouseDirs.name}, '^[0-9]')));

micronPerPixel_param1 = 0.3621;
micronPerPixel_param2 = 0.7242;

vesselTable = table();
vesselCounter = 0;

deltaPsilo = [];
deltaControl = [];
deltaMDL = [];
mouseAvgPsilo = [];
mouseAvgControl = [];
mouseAvgMDL = [];
includedMice = {};

mouseLabelsPsilo = {};
mouseLabelsControl = {};
mouseLabelsMDL = {};

allBeforeCtrl = [];
allAfterCtrl = [];
allBeforePsilo = [];
allAfterPsilo = [];
allBeforeMDL = [];
allAfterMDL = [];


%% Loop over mice
for i = 1:length(mouseDirs)
    mouseID = mouseDirs(i).name;
    mousePath = fullfile(base_folder, mouseID);

    conditionLabels = {'BeforePsilo', 'AfterPsilo', 'BeforeCtrl', 'AfterCtrl', 'BeforeMDL', 'AfterMDL'};
    folderNames = {'before-psilocybin', 'after-psilocybin', 'before-control', 'after-control', 'before-mdl-psilocybin', 'after-mdl-psilocybin'};
    treatments = {'Psilocybin', 'Psilocybin', 'Control', 'Control', 'MDL+Psilo', 'MDL+Psilo'};
    timepoints = {'Before', 'After', 'Before', 'After', 'Before', 'After'};

    allWidths = struct();

    for c = 1:length(conditionLabels)
        condPath = fullfile(mousePath, folderNames{c}, 'CapillaryVesselWidths');
        if ~isfolder(condPath)
            continue;
        end

        matFiles = dir(fullfile(condPath, '*.mat'));
        if isempty(matFiles)
            continue;
        end

        sortedFiles = sort_nat({matFiles.name});
        widths = [];
        for f = 1:length(sortedFiles)
            d = load(fullfile(condPath, sortedFiles{f}));
            if isfield(d, 'vesselWidths')
                widths = [widths; double(d.vesselWidths(:))];
            end
        end
        if isempty(widths)
            continue;
        end

        avgPix = mean(widths);
        if avgPix > 8
            scale = micronPerPixel_param1;
        else
            scale = micronPerPixel_param2;
        end
        widthsMicrons = widths * scale;
        allWidths.(conditionLabels{c}) = widthsMicrons;

        for v = 1:length(widthsMicrons)
            vesselCounter = vesselCounter + 1;
            vesselTable(vesselCounter,:) = {mouseID, treatments{c}, timepoints{c}, v, widthsMicrons(v)};
        end
    end

    if isfield(allWidths, 'BeforePsilo') && isfield(allWidths, 'AfterPsilo')
        bfPsilo = allWidths.('BeforePsilo');
        afPsilo = allWidths.('AfterPsilo');
        nP = min(length(bfPsilo), length(afPsilo));
        deltaP = 100 * (afPsilo(1:nP) - bfPsilo(1:nP)) ./ bfPsilo(1:nP);
        deltaPsilo = [deltaPsilo; deltaP];
        mouseLabelsPsilo = [mouseLabelsPsilo; repmat({mouseID}, nP, 1)];
        allBeforePsilo = [allBeforePsilo; bfPsilo(1:nP)];
        allAfterPsilo = [allAfterPsilo; afPsilo(1:nP)];
        mouseAvgPsilo(end+1) = mean(deltaP);
    end

    if isfield(allWidths, 'BeforeCtrl') && isfield(allWidths, 'AfterCtrl')
        bfCtrl = allWidths.('BeforeCtrl');
        afCtrl = allWidths.('AfterCtrl');
        nC = min(length(bfCtrl), length(afCtrl));
        deltaC = 100 * (afCtrl(1:nC) - bfCtrl(1:nC)) ./ bfCtrl(1:nC);
        deltaControl = [deltaControl; deltaC];
        mouseLabelsControl = [mouseLabelsControl; repmat({mouseID}, nC, 1)];
        allBeforeCtrl = [allBeforeCtrl; bfCtrl(1:nC)];
        allAfterCtrl = [allAfterCtrl; afCtrl(1:nC)];
        mouseAvgControl(end+1) = mean(deltaC);
    end

    if isfield(allWidths, 'BeforeMDL') && isfield(allWidths, 'AfterMDL')
        bfMDL = allWidths.('BeforeMDL');
        afMDL = allWidths.('AfterMDL');
        nM = min(length(bfMDL), length(afMDL));
        deltaM = 100 * (afMDL(1:nM) - bfMDL(1:nM)) ./ bfMDL(1:nM);
        deltaMDL = [deltaMDL; deltaM];
        mouseLabelsMDL = [mouseLabelsMDL; repmat({mouseID}, nM, 1)];
        allBeforeMDL = [allBeforeMDL; bfMDL(1:nM)];
        allAfterMDL = [allAfterMDL; afMDL(1:nM)];
        mouseAvgMDL(end+1) = mean(deltaM);
    end

    includedMice{end+1} = mouseID;
end

vesselTable.Properties.VariableNames = {'MouseID', 'Treatment', 'Timepoint', 'VesselNumber', 'VesselWidth'};


%% Plotting

% Boxplot Data of pre/post
% maxLen = max([length(allBeforeCtrl), length(allAfterCtrl), length(allBeforePsilo), length(allAfterPsilo), length(allBeforeMDL), length(allAfterMDL)]);
% padTo = @(x) [x; nan(maxLen - length(x), 1)];
% combinedVals = [padTo(allBeforeCtrl), padTo(allAfterCtrl), padTo(allBeforePsilo), padTo(allAfterPsilo), padTo(allBeforeMDL), padTo(allAfterMDL)];
% 
% subplot(1,2,1)
% hold on
% x = kron(1:6, ones(maxLen,1));
% 
% for j = 1:6
%     swarmchart(x(:,j), combinedVals(:,j), 5, 0.8*colors(j,:), 'filled', 'o', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'none');
% end
% 
% bp = boxplot(combinedVals, 'Symbol', '');
% set(bp(6,:), 'Color', 'k', 'LineWidth', 1.5);
% set(bp(5,:), 'Color', 'k', 'LineWidth', 1.5);
% xlim([0.5 6.5])
% ylabel('Vessel Width (\mum)')
% xticks(1:6)
% xticklabels({'pre','post','pre','post','pre','post'})
% t = text([1.5 3.5 5.5], repmat(min(ylim)-0.5,1,3), {'control','psilocybin','mdl+psilo'});
% set(t, 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontName', 'Arial')
% set(gca, 'FontSize', 12, 'FontName', 'Arial')
% box off

do_stats = true;
do_plots = true;

box_plot_size = [160 80];
baseline_ylimits = [-60 60];
lme_delta_ps = nan(3,3);
font_size = 4.8;
ylabels_boxplot = 'baseline % \DeltaW/W';

% Violin and Box Plots Combined for deltas
allDeltas = [deltaControl; deltaPsilo; deltaMDL];
groups = [repmat({'Control'}, numel(deltaControl), 1); repmat({'Psilocybin'}, numel(deltaPsilo), 1); repmat({'MDL+Psilo'}, numel(deltaMDL), 1)];
mouseLabels = [mouseLabelsControl; mouseLabelsPsilo; mouseLabelsMDL];

% distributionPlot({deltaControl, deltaPsilo, deltaMDL}, 'xNames', {'Control', 'Psilocybin', 'MDL+Psilo'}, 'showMM', 2, 'histOpt', 1, 'color', {colors(1,:), colors(4,:), colors(6,:)}, 'xValues', 1:3);
all_deltas_array = nan(582,3);
all_deltas_array(1:length(deltaControl),1) = deltaControl;
all_deltas_array(1:length(deltaPsilo),2) = deltaPsilo;
all_deltas_array(1:length(deltaMDL),3) = deltaMDL;

figure('Position',[100 100 box_plot_size])
plot([0 4], [0 0], '--k')
hold on
violin(all_deltas_array)

bp = boxplot(allDeltas, groups,'Symbol','');
% boxplot(allDeltas, groups, 'Colors', 'k', 'Symbol', '', 'Widths', 0.5);
set(bp(6,:), 'Color', 'k'); %median
set(bp(5,:), 'Color', 'k'); %box
ylim(baseline_ylimits)
xlim([0 4])
set(gca,'FontSize',font_size)
set(gca,'FontName','Arial')
ylabel(ylabels_boxplot,'FontSize',font_size)
xticks([1 2 3])
xticklabels({'ctrl','psil','psil+MDL'})
ax = gca;
ax.Clipping = "off";
box off

print(gcf,'-vector','-dsvg',fullfile(base_folder,'2p baseline capillaries.svg'))


%%
%run LME for this dataset
datatable = vesselTable;
datatable.Treatment = categorical(datatable.Treatment);
datatable.VesselNumber = categorical(datatable.VesselNumber);
datatable.MouseID = categorical(datatable.MouseID);
datatable.Timepoint = categorical(datatable.Timepoint);
datatable.Treatment = reordercats(datatable.Treatment, {'Psilocybin','Control', 'MDL+Psilo'});
datatable.Timepoint = reordercats(datatable.Timepoint, {'Before', 'After'});
datatable_r = datatable;
lme_model = fitlme(datatable_r, 'VesselWidth ~ Treatment*Timepoint + (1|MouseID) + (1|VesselNumber:MouseID)');
% ANOVA p-value for treatment effect
% anova_results = anova(lme_model); 
% lme_delta_ps(cur_comparison,:) = lme_model.Coefficients.pValue(2);

% disp([num2str(dt) 'AIC:']); 
% disp(lme_model.ModelCriterion.AIC);

%%

%inspect lme model resiudals
lme_residuals = residuals(lme_model);
disp(skewness(residuals(lme_model)));
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




%% Display info
disp('Preview of compiled vessel width table:');
disp(vesselTable(1:min(10,height(vesselTable)),:));

fprintf('Number of included mice: %d\n', numel(includedMice));
fprintf('Number of vessels (control): %d\n', length(deltaControl));
fprintf('Number of vessels (psilocybin): %d\n', length(deltaPsilo));
fprintf('Number of vessels (mdl+psilo): %d\n', length(deltaMDL));
