

%% set base directory
base_folder = 'D:\SN Lab\Psilocybin\BOLD simulation\';


%% (STEP 1) set general figure parameters

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
font_size = 4.8;

tiny_plot_size = [60 70];
vector_plot_size = [80 70];
box_plot_size = [120 85];
raster_plot_size = [130 250];
fc_plot_size = [250 250];


%% (STEP 2) load connectivity matrix and visualize in a circular diagram

%load data
conn_matrix = readmatrix(fullfile(base_folder,'gw_connectivity_matrix.csv'));

%normalize connectivity for visualization
conn_matrix = conn_matrix/max(conn_matrix(:));

%create circular layout
n_nodes = size(conn_matrix, 1);
angles = linspace(0, 2*pi, n_nodes+1);
angles = angles(1:end-1);
x = cos(angles);
y = sin(angles);

%create figure
figure('Color', 'white', 'Position', [100 100 100 100]);
hold on; axis equal; axis off;

%plot connections 
threshold = 0.05; % only show connections above normalization threshold
[src, tgt] = find(conn_matrix > threshold);
line_widths = 1 * conn_matrix(conn_matrix > threshold);
for k = 1:length(src)
    plot([x(src(k)), x(tgt(k))], [y(src(k)), y(tgt(k))], ...
        'Color', [0 0 0 line_widths(k)]);%, 'LineWidth', line_widths(k));
end

%plot nodes
scatter(x(1:40), y(1:40), 5, 'filled', 'MarkerFaceColor', SN_colors(1,:)); % black, left hemisphere
scatter(x(41:end), y(41:end), 5, 'filled', 'MarkerFaceColor', SN_colors(3,:)); % grey, right hemisphere

% print(gcf,'-vector','-dsvg','C:\Users\misaa\Desktop\FC sim circle.svg')


%% (STEP 3) load and plot example data (repeat for pre and post psilocybin)

%%%%%%%%%%%%
example_num = 2; %2 = pre-psilocybin, 4=post-psilocybin
%%%%%%%%%%%%
%example nums 1, 3, and 5 are simulations using linear
%extrapolation/interpolation of the pre (2) and post (4) NVC params


%plot all neural activity as raster plot
activity = readtable(fullfile(base_folder,'simulation data',num2str(example_num),'neural_activity_seed_1.csv'));

tmp = table2array(activity(2:end,2:end));
x = table2array(activity(2:end,1))/1000;
tmp2 = nan(60,80);
x2 = 1:60;
for i = 1:80
    tmp2(:,i) = interp1(x,tmp(:,i),x2,"linear");
end
tmp2 = tmp2'/30;
% tmp2 = kron(tmp2,ones(2,1));
tmp3 = repmat(tmp2,[1 1 3]);
for i = 1:size(tmp2,1)
    for j = 1:size(tmp2,2)
        if tmp2(i,j)>0
            colorval = tmp2(i,j);
            cur_color = (SN_colors(5,:)*colorval + brightwhite*(1-colorval)); %0=white, 1=pure color
        else
            colorval = max([0 tmp2(i,j)]);
            cur_color = (SN_colors(10,:)*colorval + brightwhite*(1-colorval)); %0=white, 1=pure color
        end
        tmp3(i,j,1) = cur_color(1);
        tmp3(i,j,2) = cur_color(2);
        tmp3(i,j,3) = cur_color(3);
    end
end

figure('Position',[100 100 130 130])
imshow(tmp3)
axis on
box off
% xticks([10 60 120])
% xlim([10 120])
xticks([5 30 60])
xlim([5 60])
xticklabels({'5','30','60'})
yticks([1 80])
ylim([1 80])
% yticks([1 160])
% ylim([1 160])
yticklabels({'1','80'})
set(gca,'FontSize',font_size)
set(gca,'FontName','Arial')
xlabel('time (s)','FontSize',font_size)
ylabel('brain area','FontSize',font_size)

% print(gcf,'-vector','-dsvg',fullfile(base_folder,['FC sim a' num2str(example_num) '.svg']))


% create another figure, with a colorbar
figure('Position',[100 100 100 100])
imshow(tmp3(1:end/2,1:end/2))
axis on
box off
map = interp1([0 1],[1 0.8431],0:0.01:1,'linear');
map = map';
map(:,2) = interp1([0 1],[1 0.1882],0:0.01:1,'linear');
map(:,3) = interp1([0 1],[1 0.1529],0:0.01:1,'linear');
colormap(map)
colorbar(gca,"east")


%plot BOLD for a select # of brain areas
bold = readtable(fullfile(base_folder,'simulation data',num2str(example_num),'bold_signals_seed_1.csv'));
x = table2array(bold(2:end,1))/1000;
regions_to_plot = [1 11 21 31 41];
region_colors = SN_colors([6 7 8 9 10],:);

figure('Position',[100 100 50 70])
plot(x,zeros(size(x)),'--k','Color',SN_colors(3,:))
hold on
for i = 1:length(regions_to_plot)
    tmp = table2array(bold(2:end,regions_to_plot(i)+1));
    plot(x,tmp,'-k','Color',region_colors(i,:),'LineWidth',1)
end

xlim([5 60])
ylim([-0.01 0.07])
xticks([5 30 60])
yticks([0 0.05])
box off
set(gca,'FontSize',font_size)
set(gca,'FontName','Arial')
xlabel('time (s)','FontSize',font_size)
ylabel('BOLD','FontSize',font_size)

% print(gcf,'-vector','-dsvg',fullfile(base_folder,['FC sim b' num2str(example_num) '.svg']))


%plot FC matrices
fc = readtable(fullfile(base_folder,'simulation data',num2str(example_num),'fc_results.csv'));
% fc_data = table2array(fc(1,3:end));
fc_data = mean(table2array(fc(1:20,3:end)));
fc_data = 2*(fc_data - 0.5); %0 to 1 becomes -1 to 1
% fc_data(fc_data<0) = 0;
fc_data = reshape(fc_data,[80 80]);
fc_data = kron(fc_data,ones(2,2));
fc_rgb = nan([size(fc_data) 3]);
for i = 1:size(fc_data,1)
    for j = 1:size(fc_data,2)
        if fc_data(i,j)>0
            colorval = fc_data(i,j);
            cur_color = (SN_colors(5,:)*colorval + brightwhite*(1-colorval)); %0=white, 1=pure color
        else
            colorval = -fc_data(i,j);
            cur_color = (SN_colors(10,:)*colorval + brightwhite*(1-colorval)); %0=white, 1=pure color
        end
        fc_rgb(i,j,1) = cur_color(1);
        fc_rgb(i,j,2) = cur_color(2);
        fc_rgb(i,j,3) = cur_color(3);
    end
end

figure('Position',[100 100 fc_plot_size])
imshow(fc_rgb)

xlim([1 size(fc_data,2)])
ylim([1 size(fc_data,1)])
yticks([1 size(fc_data,1)])
xticks([1 size(fc_data,2)])
% xticklabels{'1',num2str(size(fc_data,1))};
box off
axis on
set(gca,'Clipping','off')
set(gca,'YDir','reverse')
set(gca,'FontSize',font_size)
set(gca,'FontName','Arial')
xlabel('region 1','FontSize',font_size)
ylabel('region 2','FontSize',font_size)

% print(gcf,'-vector','-dsvg',fullfile(base_folder,['FC sim c' num2str(example_num) '.svg']))



%% (STEP 4) plot relationship of mean FC pairs

figure('Position',[100 100 box_plot_size])
plot([0 1],[0 1],'--k')
hold on

%get pre condition
fc = readtable(fullfile(base_folder,'simulation data',num2str(2),'fc_results.csv'));
fc_data_2 = mean(table2array(fc(1:20,3:end)));
fc_data_2 = reshape(fc_data_2,[80 80]);
fc_data_2 = triu(fc_data_2,1);
fc_data_2(fc_data_2==0) = nan;

%get post condition
fc = readtable(fullfile(base_folder,'simulation data',num2str(4),'fc_results.csv'));
fc_data_4 = mean(table2array(fc(1:20,3:end)));
fc_data_4 = reshape(fc_data_4,[80 80]);
fc_data_4 = triu(fc_data_4,1);
fc_data_4(fc_data_4==0) = nan;

scatter(fc_data_2(:),fc_data_4(:),2,'filled','MarkerFaceColor',SN_colors(10,:),'MarkerFaceAlpha',1)
xlim([0 1])
ylim([0 1])
set(gca,'FontSize',font_size)
set(gca,'FontName','Arial')
xlabel('pre-psilocybin correlation','FontSize',font_size)
ylabel('post-psilocybin correlation','FontSize',font_size)
box off

% print(gcf,'-vector','-dsvg',fullfile(base_folder,'FC sim d.svg'))



%% (STEP 5) create violin plots
extended_param_color = 0.6*ones(1,3);
series_colors = [extended_param_color;...
    SN_colors(9,:);...
    extended_param_color;...
    SN_colors(10,:);...
    extended_param_color];

figure('Position',[100 100 150 85])
% plot([0 6],[0 0],'--k')
% hold on
% x = repmat([1 2 3 4 5],[20 5]);
all_y = nan(80*80,5);
for i = 1:5 %only 2 and 4 (pre and post) are shown in the manuscript
    fc = readtable(fullfile(base_folder,'simulation data',num2str(i),'fc_results.csv'));
    fc_data = mean(table2array(fc(1:20,3:end)));
    fc_data = reshape(fc_data,[80 80]);
    fc_data = triu(fc_data,1);
    fc_data(fc_data==0) = nan;
    all_y(:,i) = fc_data(:);
end
violin(all_y)
% set(bp(6,:), 'Color', 'k'); %median
% set(bp(5,:), 'Color', 'k'); %box
% ylim([0 1])
xlim([0 6])
set(gca,'FontSize',font_size)
set(gca,'FontName','Arial')
xticklabels({'','pre','','post',''})
ylabel('correlation','FontSize',font_size)
xlabel('','FontSize',font_size)
ax = gca;
ax.Clipping = "off";
box off

% print(gcf,'-vector','-dsvg',fullfile(base_folder,'FC sim e.svg'))


%% create boxplots
extended_param_color = 0.6*ones(1,3);
series_colors = [extended_param_color;...
    SN_colors(9,:);...
    extended_param_color;...
    SN_colors(10,:);...
    extended_param_color];

figure('Position',[100 100 box_plot_size])
plot([0 6],[0 0],'--k')
hold on
x = repmat([1 2 3 4 5],[20 5]);
all_y = nan(20,5);
for i = 1:5
    fc = readtable(fullfile(base_folder,'simulation data',num2str(i),'fc_results.csv'));
    all_y(:,i) = table2array(fc(1:20,2));
    swarmchart(x(:,i),all_y(:,i),1,0.8*series_colors(i,:),'filled','o','MarkerFaceAlpha',1,'MarkerEdgeColor','none')
end
bp = boxplot(all_y,'Symbol','');
set(bp(6,:), 'Color', 'k'); %median
set(bp(5,:), 'Color', 'k'); %box
ylim([0.35 0.85])
xlim([0 6])
set(gca,'FontSize',font_size)
set(gca,'FontName','Arial')
% xticklabels({'','pre','','post',''})
ylabel('','FontSize',font_size)
xlabel('','FontSize',font_size)
ax = gca;
ax.Clipping = "off";
box off

% print(gcf,'-vector','-dsvg',fullfile(base_folder,'FC sim e.svg'))


%% (STEP 6) create parameter matrix
m = readmatrix(fullfile(base_folder,'simulation data\fc_matrix.csv'),'Range','A2:D981');
k_matrix = reshape(m(:,2),[7 7 20]);
gamma_matrix = reshape(m(:,3),[7 7 20]);
fc_matrix = reshape(m(:,4),[7 7 20]);

k_vals = mean(mean(k_matrix,3),1);
gamma_vals = mean(mean(gamma_matrix,3),2);

data = mean(fc_matrix,3);
[x, y] = meshgrid(k_vals, gamma_vals);
[xi, yi] = meshgrid(k_vals(1):0.05:k_vals(end), gamma_vals(1):0.05:gamma_vals(end)); % finer resolution

% interpolate data
zi = interp2(x, y, data, xi, yi, 'spline');

map = interp1([0 0.5 1],[0.2706 1 0.8431],0:0.01:1,'linear');
map = map';
map(:,2) = interp1([0 0.5 1],[0.4588 1 0.1882],0:0.01:1,'linear');
map(:,3) = interp1([0 0.5 1],[0.7059 1 0.1529],0:0.01:1,'linear');
% colormap(map)

%create figure, including human NVC params
figure('Position',[100 100 120 120])
contourf(xi, yi, zi, 20, 'LineColor', 'none'); % Filled contours
colormap(map);
c = colorbar;
clim([0.5 1])
c.Ticks = [0.5 1];
hold on;
scatter(0.49, 1.35, 15, 'ok','filled','MarkerFaceColor',SN_colors(9,:),'MarkerEdgeColor','k') %pre mouse
scatter(0.84, 1, 15, 'k','diamond','MarkerFaceColor',SN_colors(10,:),'MarkerEdgeColor','k') %post mouse
scatter(0.65, 0.41, 15, 'xk') %human
axis equal
xticks([0.25 1.15])
xlim([0.25 1.15])
yticks([0.25 1.75])
ylim([0.25 1.75])
box off
set(gca,'FontSize',font_size)
set(gca,'FontName','Arial')
ylabel('K','FontSize',font_size)
xlabel('Gamma','FontSize',font_size)

% print(gcf,'-vector','-dsvg',fullfile(base_folder,'FC matrix a.svg'))


%add the legend
figure('Position',[100 100 120 120])
contourf(xi, yi, zi, 20, 'LineColor', 'none'); % Filled contours
colormap(map);
c = colorbar;
clim([0.5 0.75])
c.Ticks = [0.5 0.75];
hold on;
scatter(0.49, 1.35, 15, 'ok','filled','MarkerFaceColor',SN_colors(9,:),'MarkerEdgeColor','k') %pre mouse
scatter(0.84, 1, 15, 'k','diamond','MarkerFaceColor',SN_colors(10,:),'MarkerEdgeColor','k') %post mouse
scatter(0.65, 0.41, 15, 'xk') %human
axis equal
xticks([0.25 1.15])
xlim([0.25 1.15])
yticks([0.75 1.75])
ylim([0.75 1.75])
box off
legend({'','mouse pre','mouse post','human'},'FontSize',font_size)
set(gca,'FontSize',font_size)
set(gca,'FontName','Arial')
xlabel('K','FontSize',font_size)
ylabel('Gamma','FontSize',font_size)

% print(gcf,'-vector','-dsvg',fullfile(base_folder,'FC matrix b.svg'))


%only include mouse pre/post
figure('Position',[100 100 150 150])
contourf(xi, yi, zi, 20, 'LineColor', 'none'); % Filled contours
colormap(map);
c = colorbar;
clim([0.45 0.75])
c.Ticks = [0.45 0.75];
hold on;
scatter(0.49, 1.35, 15, 'ok','filled','MarkerFaceColor',SN_colors(9,:),'MarkerEdgeColor','k') %pre mouse
scatter(0.84, 1, 15, 'k','diamond','MarkerFaceColor',SN_colors(10,:),'MarkerEdgeColor','k') %post mouse
scatter(0.65, 0.41, 15, 'xk') %human
axis equal
xticks([0.25 1.15])
xlim([0.25 1.15])
yticks([0.75 1.75])
ylim([0.75 1.75])
box off
set(gca,'FontSize',font_size)
set(gca,'FontName','Arial')
ylabel('K','FontSize',font_size)
xlabel('Gamma','FontSize',font_size)

% print(gcf,'-vector','-dsvg',fullfile(base_folder,'FC matrix c.svg'))

