% Authors: Aleah Sommers, Neosha Narayanan, January 2025
% Generates mesh, loads in parameter file, sets hydrology parameters, runs
% winter spin up equilibration for specified amount of time

clear all


steps=[1:3];

disp('Running standalone winter spinup of SHAKTI model')

if any(steps==1)
	disp('	Step 1: Mesh');

	%Generate unstructured mesh
    md = triangle(model, '../Exp/neosha17.exp', 40);

    % Optional save
	%save /jumbo/ice/ISSM/issmaleah/trunk-jpl/neosha/Models/StandaloneHydrology/neosha_Mesh17 md
end

if any(steps==2)
	disp('	Step 2: Parameterization');
	%md=loadmodel('/jumbo/ice/ISSM/issmaleah/trunk-jpl/neosha/Models/StandaloneHydrology/neosha_Mesh17');

	md=setmask(md,'','');


	% Run parameterization script to set up geometry, velocity, material properties, etc.
	md=parameterize(md,'../Par/neosha_standalone.par');

	% HYDROLOGY SPECIFIC PARAMETERIZATION:
	% Change hydrology class to Sommers' SHAKTI model
	md.hydrology=hydrologyshakti();
	
	% Define initial water head such that water pressure is 50% of ice overburden pressure
	md.hydrology.head = 0.5*md.materials.rho_ice/md.materials.rho_freshwater*md.geometry.thickness + md.geometry.base;

	% Initial subglacial gap height of 0.01m everywhere
	md.hydrology.gap_height = 0.001*ones(md.mesh.numberofelements,1);

	% Typical bed bump spacing (2m)
	md.hydrology.bump_spacing = 2.0*ones(md.mesh.numberofelements,1);

	% Typical bed bump height (0.1m)
	md.hydrology.bump_height = 0.0*ones(md.mesh.numberofelements,1);

	% Define distributed englacial input to the subglacial system (m/yr)
	% Change the value 0.0 to add distributed input
	md.hydrology.englacial_input = 0.0*ones(md.mesh.numberofvertices,1);

	% Initial Reynolds number (start at Re=1000 everywhere)
	md.hydrology.reynolds= 1000*ones(md.mesh.numberofelements,1);
    
    % Glacier front b.c.
    md.hydrology.spchead = NaN(md.mesh.numberofvertices,1);
    pos=find(md.mesh.vertexonboundary & md.mesh.y<=4.022475e6); 
    md.hydrology.spchead(pos)=md.geometry.base(pos); % Set atmospheric pressure b.c. 
    
    % Optional save
	% save /jumbo/ice/ISSM/issmaleah/trunk-jpl/neosha/Models/StandaloneHydrology/neosha_Param17 md;
end
if any(steps==3);
	disp('	Step 3: Hydro');
	%md=loadmodel('/jumbo/ice/ISSM/issmaleah/trunk-jpl/neosha/Models/StandaloneHydrology/neosha_Param17');

	md.transient=deactivateall(md.transient);
	md.transient.ishydrology=1;

	% Specify that you want to run the model on your current computer
	% Change the number of processors according to your machine
	md.cluster=generic('np', 4);

	% Define the time stepping scheme
	md.timestepping.time_step=21600/md.constants.yts; % Time step (in years), 1800=1800s=30min
	md.timestepping.final_time=730/365; % Final time (in years), 30/365=30 days
	md.settings.output_frequency=4;

    % Add moulin inputs
    time=0:md.timestepping.time_step:md.timestepping.final_time;
	md.hydrology.moulin_input = zeros(md.mesh.numberofvertices+1,numel(time));
	md.hydrology.moulin_input(end,:)=time;
    % *** this needs to be added to code for project to use random moulins ***
%     moul=find(md.mesh.vertexonboundary==0 & md.geometry.surface<2800); % Find vertices < 2900 m surface elevation
%     zi=randi(length(moul),1,10); % Select 10 random vertices for moulins
%     pos=moul(zi);
%     md.hydrology.moulin_input(pos,:)=1; % Input rate(m3/s) into moulins (steady)
    % *** end new code section for project
    
	% Specify no-flux Type 2 boundary conditions on all edges (except
	% the Type 1 condition set at the outflow above)
	md.hydrology.neumannflux=zeros(md.mesh.numberofelements+1,numel(time));
	md.hydrology.neumannflux(end,:)=time;
    
    % *** Requested outputs (melt components: meltrate, frictional, dissipation, PMPheat) ***
    %md.transient.requested_outputs={'Dummy' 'EsaEmotion' 'EsaNmotion' 'EsaUmotion'};
    % ***
    
    md.stressbalance.maxiter=100;
    
    % *** TEST convergence tolerance ***
    md.stressbalance.restol=0.05;
    md.stressbalance.reltol=0.05;
    md.stressbalance.abstol=NaN;
    md.settings.solver_residue_threshold=NaN;

	md.verbose.solution=1;
	md=solve(md,'Transient');
    
	save('../Models/winter_spin_up_17.mat', 'md', '-v7.3')

    f=md.materials.rho_freshwater./(md.materials.rho_ice.*md.geometry.thickness).*(md.results.TransientSolution(end).HydrologyHead-md.geometry.base); % Fraction of overburden
    Re=abs(md.results.TransientSolution(end).HydrologyBasalFlux)./1.787e-6; % Reynolds number

end
