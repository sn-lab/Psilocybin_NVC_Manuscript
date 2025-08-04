

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
vector_plot_size = [120 85];
box_plot_size = [120 85];
raster_plot_size = [130 250];


%% (STEP 2) plot example simulated activity, flow, BOLD
%set time vector
fps = 10;
pretime = 3;
time_vector = ((1:200)/fps)-pretime;
xlims = [-2 13];
%simulate neural activity (square wave from 0-6 s)
dt = 0.1;
x = time_vector';
Z = zeros(size(time_vector));
Z(time_vector>=0 & time_vector<=6) = 0.045;
N = size(Z, 1);

% default Balloon-Windkessel model parameters
rho = 0.34;    % Capillary resting net oxygen extraction
alpha = 0.32;  % Grubb's vessel stiffness exponent
V0 = 0.02;     % Resting blood volume fraction
k1 = 7 * rho;
k2 = 2.0;
k3 = 2 * rho - 0.2;
Gamma = 0.41*ones(N, 1);  % Rate constant for autoregulatory feedback by blood flow
K = 0.65*ones(N, 1);      % Vasodilatory signal decay
Tau = 0.98 * ones(N, 1);    % Transit time
z = 0.062; %input activity strength

% Initialize variables
X = zeros(N, 1);  % Vasodilatory signal
F = ones(N, 1);  % Blood flow
Q = ones(N, 1);  % Deoxyhemoglobin
V = ones(N, 1);  % Blood volume
BOLD = zeros(size(Z));

[BOLD, Xfull, Ffull, Qfull, Vfull] = integrateBOLD(BOLD, X, Q, F, V, Z, dt, N, rho, alpha, V0, k1, k2, k3, Gamma, K, Tau);

line_width = 1;

%plot neural activity
figure('Position',[100 100 tiny_plot_size])
tmp = Z;
tmp = tmp/max(tmp);
plot(time_vector,tmp,'k','LineWidth',line_width)
axis off
xlim(xlims)
ylim([-0.2 1.2])

% print(gcf,'-vector','-dsvg',fullfile(base_folder,'BOLD ex 1.svg'))


%plot flow
figure('Position',[100 100 tiny_plot_size])
tmp = Z;
tmp = tmp/max(tmp);
plot(time_vector,tmp,'k','Color',[0.7 0.7 0.7],'LineWidth',line_width)
hold on
tmp = Ffull-1;
tmp = tmp/max(tmp);
plot(time_vector,tmp,'k','LineWidth',line_width)
axis off
xlim(xlims)
ylim([-0.2 1.2])

% print(gcf,'-vector','-dsvg',fullfile(base_folder,'BOLD ex 2.svg'))


%plot BOLD
figure('Position',[100 100 tiny_plot_size])
tmp = Z;
tmp = tmp/max(tmp);
plot(time_vector,tmp,'k','Color',[0.7 0.7 0.7],'LineWidth',line_width)
hold on
tmp = Ffull-1;
tmp = tmp/max(tmp);
plot(time_vector,tmp,'k','Color',[0.7 0.7 0.7],'LineWidth',line_width)
tmp = BOLD;
tmp = tmp/max(tmp);
plot(time_vector,tmp,'k','LineWidth',line_width)
axis off
xlim(xlims)
ylim([-0.2 1.2])

% print(gcf,'-vector','-dsvg',fullfile(base_folder,'BOLD ex 3.svg'))


%% (STEP 3) load average vessel flow velocity change from 2P linescans and widefield ICT

%load widefield data
load(fullfile(base_folder,'speed_avgs_widefield.mat'))

%speed_means: [pts, roitype(4), cond(4)]
prepsilo_artery_avg = speed_means(:,1,3)'; %average of arteries, pre-psilocybin
postpsilo_artery_avg = speed_means(:,1,4)'; %average of arteries, post-psilocybin
prepsilo_vein_avg = speed_means(:,2,3)'; %average of arteries, pre-psilocybin
postpsilo_vein_avg = speed_means(:,2,4)'; %average of arteries, post-psilocybin

%convert from % delta to fractional change
prepsilo_artery_avg = (prepsilo_artery_avg/100)+1; 
postpsilo_artery_avg = (postpsilo_artery_avg/100)+1;
prepsilo_vein_avg = (prepsilo_vein_avg/100)+1;
postpsilo_vein_avg = (postpsilo_vein_avg/100)+1;

%set 2p time vector
P_line = 0.006; %linescan period (s)
P = P_line*50; %period between sets of 50 linescans (size of blocks for radon transformation)
t_2p = -10*P:P:66*P; %~ -3 s to +19.8 s
xlims = [-3 15];

%load 2P data
load(fullfile(base_folder,'vel_avgs_2p.mat'))

%vel_means: [pts, cond(4)]
prepsilo_capillary_avg = vel_means(:,3);
postpsilo_capillary_avg = vel_means(:,4);

%convert to fractional change of baseline (from % change)
prepsilo_capillary_avg = (prepsilo_capillary_avg/100)+1;
postpsilo_capillary_avg = (postpsilo_capillary_avg/100)+1;

%upsample to widefield time vector
prepsilo_capillary_avg = interp1(t_2p,prepsilo_capillary_avg,time_vector,"linear");
postpsilo_capillary_avg = interp1(t_2p,postpsilo_capillary_avg,time_vector,"linear");

%plot individual vessels
figure('Position',[100 100 vector_plot_size])
patch([0 6 6 0],[0 0 2.5 2.5],'k','FaceAlpha',0.15,'EdgeColor','none')
text(3,1.15,'stimulus','Color',[0.2 0.2 0.2],'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',font_size)
hold on
plot(time_vector,ones(size(time_vector)),'--k','Color',SN_colors(3,:))
plot(time_vector,prepsilo_artery_avg,'--k','Color',SN_colors(5,:))
plot(time_vector,postpsilo_artery_avg,'k','Color',SN_colors(5,:))
plot(time_vector,prepsilo_vein_avg,'--k','Color',SN_colors(10,:))
plot(time_vector,postpsilo_vein_avg,'k','Color',SN_colors(10,:))
plot(time_vector,prepsilo_capillary_avg,'--k')
plot(time_vector,postpsilo_capillary_avg,'k')
box off
xlim(xlims)
ylim([0.93 1.15])
set(gca,'FontSize',font_size)
set(gca,'FontName','Arial')
xlabel('time (s)','FontSize',font_size)
ylabel('vessel \DeltaV/V','FontSize',font_size)

% print(gcf,'-vector','-dsvg',fullfile(base_folder,'BOLD sim 1.svg'))


%create blood-volume-weigthed vessel average
capillary_fraction = 0.56;
vein_fraction = 0.22;
artery_fraction = 0.22;
prepsilo_avg = capillary_fraction*prepsilo_capillary_avg + artery_fraction*prepsilo_artery_avg + vein_fraction*prepsilo_vein_avg;
postpsilo_avg = capillary_fraction*postpsilo_capillary_avg + artery_fraction*postpsilo_artery_avg + vein_fraction*postpsilo_vein_avg;

figure('Position',[300 100 vector_plot_size])
patch([0 6 6 0],[0 0 2.5 2.5],'k','FaceAlpha',0.15,'EdgeColor','none')
text(3,1.15,'stimulus','Color',[0.2 0.2 0.2],'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',font_size)
hold on
plot(time_vector,ones(size(time_vector)),'--k','Color',SN_colors(3,:))
plot(time_vector,prepsilo_avg,'--k','Color',SN_colors(9,:))
plot(time_vector,postpsilo_avg,'k','Color',SN_colors(10,:))
box off
xlim(xlims)
ylim([0.93 1.15])
set(gca,'FontSize',font_size)
set(gca,'FontName','Arial')
xlabel('time (s)','FontSize',font_size)
ylabel('all vessel weighted avg \DeltaV/V','FontSize',font_size)

% print(gcf,'-vector','-dsvg',fullfile(base_folder,'BOLD sim 2.svg'))



%% (STEP 4) fit balloon model parameters to estimated avg flow

%order of renamed parameters:
% rho, alpha, V0, k1, k2, k3, Gamma, K, Tau
% a,   b,     c,  d,  e,  f,  g,     h, k


%fit without z adjust
startpoints = [rho, alpha, V0, k1, k2, k3, Gamma, K, Tau];
fitopts = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',startpoints*0.1,...
               'Upper',startpoints*100,...
               'StartPoint',startpoints);
ft = fittype('fitBOLD(x, a, b, c, d, e, f, g, h, k)');

%fit pre-psilocybin
y = prepsilo_avg';
[fo, gof_pre] = fit(x,y,ft,fitopts);
fo_pre = fo;
fit_flow_prepsilo = fitBOLD(x, fo.a, fo.b, fo.c, fo.d, fo.e, fo.f, fo.g, fo.h, fo.k);
params_prepsilo = [fo.g, fo.h];

%fit post-psilocybin
y = postpsilo_avg';
[fo, gof_post] = fit(x,y,ft,fitopts);
fo_post = fo;
fit_flow_postpsilo = fitBOLD(x, fo.a, fo.b, fo.c, fo.d, fo.e, fo.f, fo.g, fo.h, fo.k);
params_postpsilo = [fo.g, fo.h];


%plot Flow
figure('Position',[100 100 vector_plot_size])
patch([0 6 6 0],[0 0 1.5 1.5],'k','FaceAlpha',0.15,'EdgeColor','none')
text(3,1.15,'stimulus','Color',[0.2 0.2 0.2],'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',font_size)
hold on
plot(x,ones(size(x)),'--k','Color',SN_colors(3,:))
% plot(x,prepsilo_avg','--k','Color',colors(3,:))
plot(x,fit_flow_prepsilo,'--k','Color',colors(3,:))
% plot(x,postpsilo_avg','--k','Color',colors(4,:))
plot(x,fit_flow_postpsilo,'k','Color',colors(4,:))
xlim(xlims)
ylim([0.93 1.15])
box off
set(gca,'FontSize',font_size)
set(gca,'FontName','Arial')
xlabel('time (s)','FontSize',font_size)
ylabel('simulated \DeltaV/V','FontSize',font_size)

% print(gcf,'-vector','-dsvg',fullfile(base_folder,'BOLD sim 2.svg'))


%% (STEP 5) simulate BOLD with pre/post parameters
%run simulation with pre/post/additional parameters
Gamma_series = [nan params_prepsilo(1) nan params_postpsilo(1) nan];
K_series = [nan params_prepsilo(2) nan params_postpsilo(2) nan];
Gamma_series = fillmissing(Gamma_series,"linear");
K_series = fillmissing(K_series,"linear");

extended_param_color = 0.6*ones(1,3);
series_colors = [extended_param_color;...
    SN_colors(9,:);...
    extended_param_color;...
    SN_colors(10,:);...
    extended_param_color];

series_linewidths = [0.5 1.5 0.5 1.5 0.5];
series_linestyles = {'-','--','-','-','-'};

%plot BOLD
figure('Position',[100 100 vector_plot_size])
patch([0 6 6 0],[-1 -1 1.5 1.5],'k','FaceAlpha',0.15,'EdgeColor','none')
hold on
text(3,1.2,'stimulus','Color',[0.2 0.2 0.2],'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',font_size)
plot(x,zeros(size(x)),'--k','Color',SN_colors(3,:))

for i = [2 4] %2=pre, 4=post
    [BOLD, Xfull, Ffull, Qfull, Vfull] = integrateBOLD(BOLD, X, Q, F, V, Z, dt, N, rho, alpha, V0, k1, k2, k3, Gamma_series(i), K_series(i), Tau);
    % BOLD = BOLD/0.0029; %normalize to max of all
    plot(x,BOLD','-k','Color',series_colors(i,:),'LineStyle',series_linestyles{i})
end

xlim(xlims)
ylim([-0.001 0.004])
% ylim([-0.3 1.2])
box off
set(gca,'FontSize',font_size)
set(gca,'FontName','Arial')
xlabel('time (s)','FontSize',font_size)
ylabel('simulated BOLD','FontSize',font_size)

% print(gcf,'-vector','-dsvg',fullfile(base_folder,'BOLD sim 4.svg'))