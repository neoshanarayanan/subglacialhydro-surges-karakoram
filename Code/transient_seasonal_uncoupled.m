% Authors: Aleah Sommers, Neosha Narayanan, January 2025
% Loads in parametrized model from a previous run, loads melt timeseries, 
% and solves transient hydrology with specified timesteps

clear all;close all

% Load winter spin-up (for 2017)
load('../Models/transientmelt_2017_allmelt.mat')

% Set new initial conditions (starting from end of previous run)
md.hydrology.head=md.results.TransientSolution(end-1).HydrologyHead;
md.hydrology.gap_height=md.results.TransientSolution(end-1).HydrologyGapHeight;
md.hydrology.reynolds=md.results.TransientSolution(end-1).HydrologyBasalFlux./1.787e-6;

% 8 March 2023
% This script reads in the glacier melt data separated by year that was
% created using read_process_melt.m to process the melt data, then
% interpolates onto a SHAKTI model mesh and interpolates from daily values
% to smaller time steps

% Load glacier melt data, separated by year
load('../Data/total_melt.mat') % this file comes from read_process_melt

simname = melt17;

% Time-stepping
md.cluster=generic('np',20);

md.timestepping.time_step=1800/md.constants.yts; % Time step (in years)
md.timestepping.final_time=365/365; % number of days
timevec=0:md.timestepping.time_step:md.timestepping.final_time;
% Only save model output every day to reduce file size
md.settings.output_frequency=(md.timestepping.time_step*365)^-1; 

% Make vectors of lat/lon into 2D matrices for conversion
lat2=repmat(lat,1,length(lon));lon2=repmat(lon',length(lat),1);

% Get lat/long corresponding to our x/y coordinates (with UTM Zone 43)
[md.mesh.lat,md.mesh.long]=utm2ll(md.mesh.x,md.mesh.y,43);


% Interpolate daily melt onto mesh (m/yr) -- for a particular year
englacial_input_daily=zeros(length(md.mesh.lat),size(simname,3));
for i=1:size(simname,3)
    englacial_input_daily(:,i)=InterpFromGrid(lat2(:,1),lon2(1,:),simname(:,:,i),md.mesh.lat,md.mesh.long);
end

% Now interpolate from daily to smaller time steps
tm=size(englacial_input_daily,2);
ieb=zeros(md.mesh.numberofvertices,length(1:md.timestepping.time_step*365:tm+1));
for i=1:md.mesh.numberofvertices
    ieb(i,:)=interp1(englacial_input_daily(i,:),1:md.timestepping.time_step*365:tm+1);
end

ieb(isnan(ieb))=0; % ensure no NaNs

% Now put it into the model!
md.hydrology.englacial_input=zeros(md.mesh.numberofvertices+1,length(timevec));
md.hydrology.englacial_input(end,:)=timevec; % This last row will be a time index to keep track
md.hydrology.englacial_input(1:end-1,:)=ieb; % Put in our 30-minute melt values

% Specify 20 processors
md.cluster=generic('np',20);

% Solve the model 
md = solve(md, 'Transient');
txt = '2017, starting from end of winter spinup';
save('../Models/transientmelt_2017_allmelt.mat', 'md', 'txt')

