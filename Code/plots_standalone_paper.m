% Author: Neosha Narayanan, January 2025
% Generate plots for hydrology/surge paper

% Load data
load ../Models/StandaloneHydrology/trianglemesh_latemay/transientmelt_2017_allmelt_2.mat

%% Geometry plots

plotmodel(md, 'data', md.geometry.surface, 'unit#all','km')
colormap(brewermap([], "OrRd"))
title('Surface Elevation', '(m asl)')
ylabel('y (km)')
xlabel('x (km)')

plotmodel(md, 'data', md.geometry.base, 'unit#all','km')
colormap(brewermap([], "OrRd"))
title('Bed Elevation', '(m asl)')
ylabel('y (km)')
xlabel('x (km)')

plotmodel(md, 'data', md.geometry.thickness, 'unit#all','km')
colormap(brewermap([], 'YlOrRd'))
title('Ice Thickness', '(m)')
ylabel('y (km)')
xlabel('x (km)')

%% Meltwater inputs - plot 
clear

load ../Data/has1glaciermelt.mat
glaciermelt = melt;

load ../Data/has1snowmelt.mat
snowmelt = total_melt;

load ../Data/has1liquidprecip.mat
liquidprecip = melt;

[i1,j1]=find(glaciermelt(:,:,180)~=0);
glaciermelt_test=glaciermelt(round(mean(i1(:))),round(mean(j1(:))),:);
%glaciermelt_test = melt(i(1), j(1), :);
glaciermelt_test=squeeze(glaciermelt_test);

[i2,j2]=find(snowmelt(:,:,180)~=0);
snowmelt_test=snowmelt(round(mean(i2(:))),round(mean(j2(:))),:);
snowmelt_test=squeeze(snowmelt_test);

[i3,j3]=find(liquidprecip(:,:,180)~=0);
liquidprecip_test=liquidprecip(round(mean(i3(:))),round(mean(j3(:))),:);
liquidprecip_test=squeeze(liquidprecip_test);

figure;plot(glaciermelt_test/365*1000, 'LineWidth', 2) % This plots all five years 2017-2021 for our sample point
hold on
plot(snowmelt_test/365*1000, 'LineWidth', 2)
plot(liquidprecip_test/365*1000, 'LineWidth', 2)
ylabel('Englacial Input (mm d^{-1})')
xlabel('Date')
title('Surface Water Inputs')
xticks(1:182:1827) % Set the ticks every 6 months
dateaxis('x', 10, datetime(2016,12,31))
ax=gca;
ax.FontSize = 15;
ax.TitleFontSizeMultiplier = 1.5;
labels=string(ax.XAxis.TickLabels);
labels(1:2:end)=nan;
%labels(end) = nan;
ax.XAxis.TickLabels = labels;
legend('Ice Melt','Snow Melt','Liquid Precipitation', 'location', 'northwest')
hold off

%% Now plot melt inputs in space 

allmelt = snowmelt + glaciermelt + liquidprecip;
avgmelt = mean(allmelt(:, :, 135:304), 3);

figure
contourf(lat, lon, allmelt(:, :, 180), 50);
%title('Average Melt Season Surface Melt and Liquid Precipitation')
xlabel('Longitude (º)')
ylabel('Latitude (º)')
%imagesc(lat, lon, allmelt(:, :, 180))
view([90 -90])
colormap(brewermap([], "Blues"))
colorbar
xlim([36.34 36.465])
ylim([74.565 74.68])

a=colorbar;
a.Label.String = 'Englacial Input (mm d^{-1})';

ax=gca;
ax.FontSize = 15;
ax.TitleFontSizeMultiplier = 1.2;


%% Temporal plot of 2017

load ../Models/StandaloneHydrology/trianglemesh_latemay/transientmelt_2017_allmelt_2.mat

for tt=1:length(md.results.TransientSolution)
    bmean(tt)=mean(md.results.TransientSolution(tt).HydrologyGapHeight);
    bmin(tt)=min(md.results.TransientSolution(tt).HydrologyGapHeight);
    bmax(tt)=max(md.results.TransientSolution(tt).HydrologyGapHeight);
    
    hmean(tt)=mean(md.results.TransientSolution(tt).HydrologyHead);
    hmin(tt)=min(md.results.TransientSolution(tt).HydrologyHead);
    hmax(tt)=max(md.results.TransientSolution(tt).HydrologyHead);
    
    qmean(tt)=mean(md.results.TransientSolution(tt).HydrologyBasalFlux);
    qmin(tt)=min(md.results.TransientSolution(tt).HydrologyBasalFlux);
    qmax(tt)=max(md.results.TransientSolution(tt).HydrologyBasalFlux);
    
    Nmean(tt)=mean(md.results.TransientSolution(tt).EffectivePressure);
    Nmin(tt)=min(md.results.TransientSolution(tt).EffectivePressure);
    Nmax(tt)=max(md.results.TransientSolution(tt).EffectivePressure);
    
    f=md.materials.rho_freshwater./(md.materials.rho_ice.*md.geometry.thickness).*(md.results.TransientSolution(tt).HydrologyHead-md.geometry.base); % Fraction of overburden
    f(f<0)=0;
    fmean(tt)=mean(f);
    fmin(tt)=min(f);
    fmax(tt)=max(f);
    
end

linewidth = 1.5; 

subplot(3,1,1)
plot(1:tt,hmean,1:tt,hmin,1:tt,hmax, 'LineWidth',linewidth)
%plot(1:tt,hmean,1:tt,hmin, 'LineWidth',linewidth)
%xlim([0 366])
%xticks(1:31:365)
dateaxis('x', 3, datetime(2017,1,1))
ylabel('Head (m)')


subplot(3,1,2)
plot(1:tt,qmean,1:tt,qmin,1:tt,qmax, 'LineWidth',linewidth)
%plot(1:tt,qmean,1:tt,qmin, 'LineWidth',linewidth)
xline(213)
%xline(1)
%xline(278)
%xline(291)
xlim([0 366])
xticks(1:31:365)
%dateaxis('x', 3, datetime(2017,1,1))
ylabel('Basal flux (m^2 s^{-1})')

subplot(3,1,3)
plot(1:tt,Nmean,1:tt,Nmin,1:tt,Nmax,'LineWidth',linewidth)
%plot(1:tt,Nmean,1:tt,Nmin,'LineWidth',linewidth)
xlim([0 366])
xticks(1:31:365)
dateaxis('x', 3, datetime(2017,1,1))
ylabel({'Effective pressure'; '(Pa)'})
xlabel('Date')


sgtitle('2017 Model Outputs')
fontsize(gcf,scale=1.2)

legend('Mean','Min','Max')

%% Model results (temporal plots)
% load ../Models/StandaloneHydrology/trianglemesh_latemay/transientmelt_2017_allmelt_2.mat

% MAKE BIG PLOT 

% First, get daily melt inputs
englacial_inputs = md.hydrology.englacial_input;
meltinputs_2017 = mean(englacial_inputs(:, 1:48:end), 1)/365*1000;


% Then, get hydraulic heads 
for t=1:length(md.results.TransientSolution)
    hmeans_2017(t)=mean(md.results.TransientSolution(t).HydrologyHead);
end

figure
yyaxis left
plot(meltinputs_2017, 'LineWidth', 2)
ylabel('Melt + Precipitation Input (mm/d)')
yyaxis right
plot(hmeans_2017, 'LineWidth', 2)
ylabel('Hydraulic Head (m)')
title('Hydraulic Head Response to Melt Input')
xticks(1:31:365) % Set the ticks every 6 months
ax=gca;
ax.FontSize = 15;
ax.TitleFontSizeMultiplier = 1.5;
xlim([0 365])
dateaxis('x', 3, datetime(2017,1,1))
legend('Englacial Input Forcing', 'Hydraulic Head')

%% Now get effective pressures at June 5 and June 15: 

plotmodel(md, 'data', 1e-6 * md.results.TransientSolution(155).EffectivePressure, 'unit#all','km')
colormap(flipud(brewermap([], "Reds")))
%colormap(brewermap([], "Reds"))
title('June 5', 'Effective Pressure (MPa)')
ylabel('y (km)')
xlabel('x (km)')

lower_lim = min(1e-6 * md.results.TransientSolution(155).EffectivePressure);
upper_lim = max(1e-6 * md.results.TransientSolution(155).EffectivePressure);


plotmodel(md, 'data', 1e-6 * md.results.TransientSolution(160).EffectivePressure, 'unit#all','km')
colormap(flipud(brewermap([], "Reds")))
%colormap(brewermap([], "Reds"))
title('June 10', 'Effective Pressure (MPa)')
ylabel('y (km)')
xlabel('x (km)')
clim([lower_lim upper_lim])

%% Basal flux during the spike 

lower_lim = min(log10(md.results.TransientSolution(195).HydrologyBasalFlux));
upper_lim = max(log10(md.results.TransientSolution(195).HydrologyBasalFlux));

% plotmodel(md, 'data', log10(md.results.TransientSolution(155).HydrologyBasalFlux), 'unit#all','km')
% colormap(brewermap([], "PuBu"))
% title('June 5', 'Log_{10} Basal Flux (m^2 s^{-1})')
% ylabel('y (km)')
% xlabel('x (km)')
% clim([lower_lim upper_lim])
% 
% plotmodel(md, 'data', log10(md.results.TransientSolution(160).HydrologyBasalFlux), 'unit#all','km')
% colormap(brewermap([], "PuBu"))
% title('June 10', 'Log_{10} Basal Flux (m^2 s^{-1})')
% ylabel('y (km)')
% xlabel('x (km)')
% clim([lower_lim upper_lim])

% plotmodel(md, 'data', log10(md.results.TransientSolution(165).HydrologyBasalFlux), 'unit#all','km')
% colormap(brewermap([], "PuBu"))
% title('June 15', 'Log_{10} Basal Flux (m^2 s^{-1})')
% ylabel('y (km)')
% xlabel('x (km)')
% clim([lower_lim upper_lim])


plotmodel(md, 'data', log10(md.results.TransientSolution(180).HydrologyBasalFlux), 'unit#all','km')
colormap(brewermap([], "PuBu"))
title('June 30', 'Log_{10} Basal Flux (m^2 s^{-1})')
ylabel('y (km)')
xlabel('x (km)')
clim([lower_lim upper_lim])

% % plotmodel(md, 'data', log10(md.results.TransientSolution(200).HydrologyBasalFlux), 'unit#all','km')
% % colormap(brewermap([], "PuBu"))
% % title('July 19', 'Log_{10} Basal Flux (m^2 s^{-1})')
% % ylabel('y (km)')
% % xlabel('x (km)')
% % clim([lower_lim upper_lim])

% plotmodel(md, 'data', log10(md.results.TransientSolution(211).HydrologyBasalFlux), 'unit#all','km')
% colormap(brewermap([], "PuBu"))
% title('July 30', 'Log_{10} Basal Flux (m^2 s^{-1})')
% ylabel('y (km)')
% xlabel('x (km)')
% clim([lower_lim upper_lim])

%% Basal flux during the year 2017


lower_lim = min(log10(md.results.TransientSolution(260).HydrologyBasalFlux));
upper_lim = max(log10(md.results.TransientSolution(260).HydrologyBasalFlux));

% plotmodel(md, 'data', log10(md.results.TransientSolution(1).HydrologyBasalFlux), 'unit#all','km')
% colormap(brewermap([], "PuBu"))
% title('January 1', 'Log_{10} Basal Flux (m^2 s^{-1})')
% ylabel('y (km)')
% xlabel('x (km)')
% clim([lower_lim upper_lim])

plotmodel(md, 'data', log10(md.results.TransientSolution(213).HydrologyBasalFlux), 'unit#all','km')
colormap(brewermap([], "PuBu"))
title('August 1', 'Log_{10} Basal Flux (m^2 s^{-1})')
ylabel('y (km)')
xlabel('x (km)')
clim([lower_lim upper_lim])


% plotmodel(md, 'data', log10(md.results.TransientSolution(230).HydrologyBasalFlux), 'unit#all','km')
% colormap(brewermap([], "PuBu"))
% title('August 18', 'Log_{10} Basal Flux (m^2 s^{-1})')
% ylabel('y (km)')
% xlabel('x (km)')
% clim([lower_lim upper_lim])


% plotmodel(md, 'data', log10(md.results.TransientSolution(255).HydrologyBasalFlux), 'unit#all','km')
% colormap(brewermap([], "OrRd"))
% title('September 12', 'Log_{10} Basal Flux (m^2 s^{-1})')
% ylabel('y (km)')
% xlabel('x (km)')
% clim([lower_lim upper_lim])

% plotmodel(md, 'data', log10(md.results.TransientSolution(278).HydrologyBasalFlux), 'unit#all','km')
% colormap(brewermap([], "PuBu"))
% title('October 5', 'Log_{10} Basal Flux (m^2 s^{-1})')
% ylabel('y (km)')
% xlabel('x (km)')
% clim([lower_lim upper_lim])



% plotmodel(md, 'data', log10(md.results.TransientSolution(291).HydrologyBasalFlux), 'unit#all','km')
% colormap(brewermap([], "PuBu"))
% title('October 18', 'Log_{10} Basal Flux (m^2 s^{-1})')
% ylabel('y (km)')
% xlabel('x (km)')
% clim([lower_lim upper_lim])

% 
% plotmodel(md, 'data', log10(md.results.TransientSolution(335).HydrologyBasalFlux), 'unit#all','km')
% colormap(brewermap([], "PuBu"))
% title('December 1', 'Log_{10} Basal Flux (m^2 s^{-1})')
% ylabel('y (km)')
% xlabel('x (km)')
% clim([lower_lim upper_lim])

%% Transmissivity before and after June spike 


june20_daynumber = 171;
sept4_daynumber = 247;
oct10_daynumber = 283;

% Define a few constants
nu=1.787e-6; % Kinematic viscosity (m2/s)
omega=0.001; % Parameter controlling transition from laminar to turbulent flow

% June 20
Re=md.results.TransientSolution(june20_daynumber).HydrologyBasalFlux./nu;
% Calculate local transmissivity
K=md.results.TransientSolution(june20_daynumber).HydrologyGapHeight.^3.*md.constants.g./(12*nu.*(1+omega.*Re));
plotmodel(md,'data',log10(K.*3600*24),'unit#all','km')
title('June 20')
ylabel('y (km)')
xlabel('x (km)')
colormap(brewermap([], "YlOrRd"))
clim([0 7])

a=colorbar;
a.Label.String = 'Log_{10} Transmissivity (m^2/d)';
june20_k = log10(K.*3600*24);

% Do the same for September 4
Re=md.results.TransientSolution(sept4_daynumber).HydrologyBasalFlux./nu;
K=md.results.TransientSolution(sept4_daynumber).HydrologyGapHeight.^3.*md.constants.g./(12*nu.*(1+omega.*Re));
plotmodel(md,'data',log10(K.*3600*24),'unit#all','km')
title('September 4')
ylabel('y (km)')
xlabel('x (km)')
colormap(brewermap([], "YlOrRd"))
clim([0 7])
a=colorbar;
a.Label.String = 'Log_{10} Transmissivity (m^2/d)';
sept4_k = log10(K.*3600*24);

figure
histogram(june20_k, 20)
hold on
histogram(sept4_k, 20)
hold off
legend('June 20', 'September 4')
% october 10
% Re=md.results.TransientSolution(oct10_daynumber).HydrologyBasalFlux./nu;
% K=md.results.TransientSolution(oct10_daynumber).HydrologyGapHeight.^3.*md.constants.g./(12*nu.*(1+omega.*Re));
% plotmodel(md,'data',log10(K.*3600*24),'unit#all','km')
% title('October 10', 'Log_{10} Basal Flux (m^2 s^{-1})')
% ylabel('y (km)')
% xlabel('x (km)')
% colormap(brewermap([], "PuRd"))


%% Figure 3: 
%clear all 

% load('../Models/StandaloneHydrology/trianglemesh_latemay/transientmelt_2017_allmelt.mat')
% md2017_1 = md;
% 
% load('../Models/StandaloneHydrology/trianglemesh_latemay/transientmelt_2017_allmelt_2.mat')
% md2017_2 = md;

lower_lim = min(log10(md2017_1.results.TransientSolution(1).HydrologyBasalFlux));
upper_lim = max(log10(md2017_1.results.TransientSolution(1).HydrologyBasalFlux));

plotmodel(md2017_1, 'data', log10(md2017_1.results.TransientSolution(1).HydrologyBasalFlux), 'unit#all','km')
colormap(brewermap([], "PuBu"))
ylabel('y (km)')
xlabel('x (km)')
clim([lower_lim+1.5 upper_lim])
xlim([461.1 467])
ylim([4022.1 4034])
a=colorbar;
a.Label.String = 'Basal Flux (m^2/s)';
title('Following Equilibration Spinup', 'Log_{10} Basal Flux (m^2 s^{-1})')

plotmodel(md2017_2, 'data', log10(md2017_2.results.TransientSolution(end).HydrologyBasalFlux), 'unit#all','km')
colormap(brewermap([], "PuBu"))
ylabel('y (km)')
xlabel('x (km)')
xlim([461.1 467])
ylim([4022.1 4034])
clim([lower_lim+1.5 upper_lim])
a=colorbar;
a.Label.String = 'Basal Flux (m^2/s)';
title('Following Equilibration Spinup + Melt Cycle', 'Log_{10} Basal Flux (m^2 s^{-1})')


difference = md2017_1.results.TransientSolution(1).HydrologyBasalFlux - md2017_2.results.TransientSolution(end).HydrologyBasalFlux;
plotmodel(md2017_1, 'data', difference, 'unit#all','km')
ylabel('y (km)')
xlabel('x (km)')
colormap(flipud(brewermap([], "OrRd")))
xlim([461.1 467])
ylim([4022.1 4034])
a=colorbar;
a.Label.String = 'Difference in Basal Flux (m^2/s)';
clim([-10*1e-4 0])
title('Difference')

%% Create big surge plot again: first load parts 

clear all

load('../Models/StandaloneHydrology/trianglemesh_latemay/transientmelt_2017_allmelt_2.mat')
md2017 = md;

load('../Models/StandaloneHydrology/trianglemesh_latemay/transientmelt_2018_allmelt.mat')
md2018 = md;

load('../Models/StandaloneHydrology/trianglemesh_latemay/transientmelt_2019_allmelt.mat')
md2019 = md;

load('../Models/StandaloneHydrology/trianglemesh_latemay/transientmelt_2020_allmelt.mat')
md2020 = md;

load('../Models/StandaloneHydrology/trianglemesh_latemay/transientmelt_2021_allmelt.mat')
md2021 = md;

%% Plot all the years together
for tt=1:length(md2017.results.TransientSolution)
    
    hmean2017(tt)=mean(md2017.results.TransientSolution(tt).HydrologyHead);
    hmin2017(tt)=min(md2017.results.TransientSolution(tt).HydrologyHead);
    hmax2017(tt)=max(md2017.results.TransientSolution(tt).HydrologyHead);
    
    qmean2017(tt)=mean(md2017.results.TransientSolution(tt).HydrologyBasalFlux);
    qmin2017(tt)=min(md2017.results.TransientSolution(tt).HydrologyBasalFlux);
    qmax2017(tt)=max(md2017.results.TransientSolution(tt).HydrologyBasalFlux);
    
    Nmean2017(tt)=mean(md2017.results.TransientSolution(tt).EffectivePressure);
    Nmin2017(tt)=min(md2017.results.TransientSolution(tt).EffectivePressure);
    Nmax2017(tt)=max(md2017.results.TransientSolution(tt).EffectivePressure);
  
end

for tt=1:length(md2018.results.TransientSolution)
    
    hmean2018(tt)=mean(md2018.results.TransientSolution(tt).HydrologyHead);
    hmin2018(tt)=min(md2018.results.TransientSolution(tt).HydrologyHead);
    hmax2018(tt)=max(md2018.results.TransientSolution(tt).HydrologyHead);
    
    qmean2018(tt)=mean(md2018.results.TransientSolution(tt).HydrologyBasalFlux);
    qmin2018(tt)=min(md2018.results.TransientSolution(tt).HydrologyBasalFlux);
    qmax2018(tt)=max(md2018.results.TransientSolution(tt).HydrologyBasalFlux);
    
    Nmean2018(tt)=mean(md2018.results.TransientSolution(tt).EffectivePressure);
    Nmin2018(tt)=min(md2018.results.TransientSolution(tt).EffectivePressure);
    Nmax2018(tt)=max(md2018.results.TransientSolution(tt).EffectivePressure);
  
end

for tt=1:length(md2019.results.TransientSolution)
    
    hmean2019(tt)=mean(md2019.results.TransientSolution(tt).HydrologyHead);
    hmin2019(tt)=min(md2019.results.TransientSolution(tt).HydrologyHead);
    hmax2019(tt)=max(md2019.results.TransientSolution(tt).HydrologyHead);
    
    qmean2019(tt)=mean(md2019.results.TransientSolution(tt).HydrologyBasalFlux);
    qmin2019(tt)=min(md2019.results.TransientSolution(tt).HydrologyBasalFlux);
    qmax2019(tt)=max(md2019.results.TransientSolution(tt).HydrologyBasalFlux);
    
    Nmean2019(tt)=mean(md2019.results.TransientSolution(tt).EffectivePressure);
    Nmin2019(tt)=min(md2019.results.TransientSolution(tt).EffectivePressure);
    Nmax2019(tt)=max(md2019.results.TransientSolution(tt).EffectivePressure);
  
end

for tt=1:length(md2020.results.TransientSolution)
    
    hmean2020(tt)=mean(md2020.results.TransientSolution(tt).HydrologyHead);
    hmin2020(tt)=min(md2020.results.TransientSolution(tt).HydrologyHead);
    hmax2020(tt)=max(md2020.results.TransientSolution(tt).HydrologyHead);
    
    qmean2020(tt)=mean(md2020.results.TransientSolution(tt).HydrologyBasalFlux);
    qmin2020(tt)=min(md2020.results.TransientSolution(tt).HydrologyBasalFlux);
    qmax2020(tt)=max(md2020.results.TransientSolution(tt).HydrologyBasalFlux);
    
    Nmean2020(tt)=mean(md2020.results.TransientSolution(tt).EffectivePressure);
    Nmin2020(tt)=min(md2020.results.TransientSolution(tt).EffectivePressure);
    Nmax2020(tt)=max(md2020.results.TransientSolution(tt).EffectivePressure);
  
end

for tt=1:length(md2021.results.TransientSolution)
    
    hmean2021(tt)=mean(md2021.results.TransientSolution(tt).HydrologyHead);
    hmin2021(tt)=min(md2021.results.TransientSolution(tt).HydrologyHead);
    hmax2021(tt)=max(md2021.results.TransientSolution(tt).HydrologyHead);
    
    qmean2021(tt)=mean(md2021.results.TransientSolution(tt).HydrologyBasalFlux);
    qmin2021(tt)=min(md2021.results.TransientSolution(tt).HydrologyBasalFlux);
    qmax2021(tt)=max(md2021.results.TransientSolution(tt).HydrologyBasalFlux);
    
    Nmean2021(tt)=mean(md2021.results.TransientSolution(tt).EffectivePressure);
    Nmin2021(tt)=min(md2021.results.TransientSolution(tt).EffectivePressure);
    Nmax2021(tt)=max(md2021.results.TransientSolution(tt).EffectivePressure);
  
end


total_time = length(md2017.results.TransientSolution)+length(md2018.results.TransientSolution) +length(md2019.results.TransientSolution)+length(md2020.results.TransientSolution)+length(md2021.results.TransientSolution);
t = datetime(2017, 1, 1:total_time);

hmean = [hmean2017 hmean2018 hmean2019 hmean2020 hmean2021];
hmin = [hmin2017 hmin2018 hmin2019 hmin2020 hmin2021];
hmax = [hmax2017 hmax2018 hmax2019 hmax2020 hmax2021];

qmean = [qmean2017 qmean2018 qmean2019 qmean2020 qmean2021];
qmin = [qmin2017 qmin2018 qmin2019 qmin2020 qmin2021];
qmax = [qmax2017 qmax2018 qmax2019 qmax2020 qmax2021];

Nmean = [Nmean2017 Nmean2018 Nmean2019 Nmean2020 Nmean2021];
Nmin = [Nmin2017 Nmin2018 Nmin2019 Nmin2020 Nmin2021];
Nmax = [Nmax2017 Nmax2018 Nmax2019 Nmax2020 Nmax2021];

figure
% subplot(2, 1, 1)
% plot(t, hmean, 'LineWidth', 1.5)
% xticks(datetime(2017, 1, 1:182:1827)) % Set the ticks every 6 months
% xtickformat('MM/yyyy')
% hold on
% plot(t, hmin, 'LineWidth', 1.5)
% plot(t, hmax, 'LineWidth', 1.5)
% ylabel('Hydraulic Head (m)')
% hold off

subplot(2, 1, 1)
plot(t, qmean, 'LineWidth', 1.5)
xticks(datetime(2017, 1, 1:182:1827)) % Set the ticks every 6 months
xtickformat('MM/yyyy')
hold on
plot(t, qmin, 'LineWidth', 1.5)
plot(t, qmax, 'LineWidth', 1.5)
ylabel('Basal Flux (m^2/s)')
hold off

subplot(2, 1, 2)
plot(t, Nmean, 'LineWidth', 1.5)
xticks(datetime(2017, 1, 1:182:1827)) % Set the ticks every 6 months
xtickformat('MM/yyyy')
hold on
plot(t, Nmin, 'LineWidth', 1.5)
plot(t, Nmax, 'LineWidth', 1.5)
ylabel('Effective Pressure (Pa)')
xlabel('Date')
legend('Mean', 'Min', 'Max')
hold off


%% Now generate the same plot but with Flavien's velocity measurements on top 


vel_data = readtable('../Data/beaud-velocities-surge.csv');
vel_dates = vel_data.Date;
velocities = vel_data.Velocity;

figure
sgtitle('Surge Phases 2017-2019 at Shishper Glacier')
%colororder(flipud('sail'))




% FIRST PLOT

subplot(2, 1, 1)

ax=gca;
%ax.YAxis(1).Color = [0.47 0.25 0.80];
plot(t, qmean, 'LineWidth', 1.5)
hold on
%plot(t, qmin, 'LineWidth', 1.5)
%plot(t, qmax, 'LineWidth', 1.5)
ylabel('Basal Flux (m^2/s)')

hold off
%stairs(vel_dates+10, velocities, 'LineWidth', 1.5)
xlim([datetime('1-Jan-2017') datetime('31-Dec-2019')])
%ylabel('Velocity (m/d)')

% line = xline(datetime('23-June-2019')); % GLOF!
% line.Color = 'r';
% line.LineWidth = 2;
% line.LabelOrientation = 'horizontal';
% lgd = legend('SHAKTI-simulated hydrology', 'Beaud et al. (2023) data');
% lgd.Location = "southeast";


ax.FontSize = 12;



% SECOND

subplot(2, 1, 2)
yyaxis left
plot(t, Nmean, 'LineWidth', 1.5)
hold on
%plot(t, Nmin, 'LineWidth', 1.5)
%plot(t, Nmax, 'LineWidth', 1.5)
ylabel('Effective Pressure (Pa)')


yyaxis right
stairs(vel_dates+10, velocities, 'LineWidth', 1.5)
xlim([datetime('1-Jan-2017') datetime('31-Dec-2019')])
ylabel('Velocity (m/d)')
ax=gca;
ax.FontSize = 12;
line = xline(datetime('23-June-2019')); % GLOF!
line.Color = 'r';
line.LineWidth = 2;
line.LabelOrientation = 'horizontal';
line.Label = "GLOF (June 23, 2019)";
hold off
lgd = legend('SHAKTI-simulated hydrology', 'Beaud et al. (2023) data');
lgd.Location = "southeast";




%% Make a plot without the velocities (for editing with Illustrator)

figure
sgtitle('Surge Phases 2017-2019 at Shishper Glacier')
colororder(flipud('sail'))

subplot(2, 1, 1);
set(gca, 'Position', [0.1, 0.35, 0.8, 0.455]);
plot(t, qmean, 'LineWidth', 1.5)
ylabel('Basal Flux (m^2/s)')
xlim([datetime('1-Jan-2017') datetime('31-Dec-2019')])

subplot(2, 1, 2);
set(gca, 'Position', [0.1, 0.1, 0.8, 0.65]);
plot(t, Nmean, 'LineWidth', 1.5)
xlim([datetime('1-Jan-2017') datetime('31-Dec-2019')])
ylabel('Effective Pressure (Pa)')

line = xline(datetime('23-June-2019')); % GLOF!
line.Color = 'r';
line.LineWidth = 2;
line.LabelOrientation = 'horizontal';
line.Label = "GLOF (June 23, 2019)";
hold off

%% Try again 

figure
sgtitle('Surge Phases 2017-2019 at Shishper Glacier')

ax1 = axes('Position', [0.1, 0.625, 0.8, 0.2]); 
ax1.FontSize = 40;
ax1.YAxis.Color = [0.7 0.7 0.7];
%plot(ax1, t, qmean, 'LineWidth', 1.5, 'Color', [0.3010 0.7450 0.9330])
patch(ax1, t, qmean, [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1.5)
ylabel(ax1, 'Basal Flux (m^2/s)', 'FontSize', 14)
xlim(ax1, [datetime('1-Jan-2017') datetime('31-Dec-2019')])
set(gca,'YAxisLocation','right')
%line = xline(ax1, datetime('23-June-2019')); % GLOF!
%line.Color = 'r';
%line.LineWidth = 2;

set(ax1,'xticklabel',[])

ax2 = axes('Position', [0.1, 0.1, 0.8, 0.45]); 
ax2.FontSize = 40;
yyaxis left
plot(ax2, t, Nmean, 'LineWidth', 1.5, 'Color', [0.4940 0.1840 0.5560])
ax2.YAxis(1).Color = [0.4940 0.1840 0.5560];
xlim(ax2, [datetime('1-Jan-2017') datetime('31-Dec-2019')])
ylabel(ax2, 'Effective Pressure (Pa)', 'FontSize', 15)

yyaxis right
ax2.YAxis(2).Color = [0 0 0];
ylabel(ax2, 'Velocity (m/d)', 'FontSize', 15)
ylim(ax2, [0 7])

% Add the GLOF event line to the second plot
line = xline(ax2, datetime('23-June-2019')); % GLOF!
line.Color = 'r';
line.LineWidth = 2;
line.LabelOrientation = 'horizontal';
line.Label = "GLOF (June 23, 2019)";
line.LabelHorizontalAlignment = 'left';


set(ax1, 'FontSize', 12);
set(ax2, 'FontSize', 12);
