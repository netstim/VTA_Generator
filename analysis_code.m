%% set parameters:
reforce = 1; % set to 1 to make sure everything is recalculated

% Params for Step 1: Calculate Pseudo VTAs
amp.mean=2.5; % average vta amplitude
amp.jit=2; % jitter of vta amplitudes
numvtas=1000; % how many vtas to generate
vta_jit=1; % jitter of spatial spread
roidef=fullfile(ea_space([],'atlases'),'DISTAL Nano (Ewert 2017)','lh','STN.nii.gz'); % area to cover by VTAs

% Params for Step 2: Get Pseudo Improvements
sweetspots={fullfile(ea_space([],'atlases'),'TOR-signPD (Boutet 2021)','lh','tremor_hotspot.nii.gz'),... % sweetspot example 1 (nifti)
    fullfile(ea_space([],'atlases'),'DBS Tractography Atlas (Middlebrooks 2020)','lh','DRTT.mat'),... % sweetspot example 2 (tract))
    fullfile(ea_space([],'atlases'),'TOR-signPD (Boutet 2021)','lh','tremor_coldspot.nii.gz'),... % sweetspot example 1 (nifti)
    fullfile(ea_space([],'atlases'),'DBS Tractography Atlas (Middlebrooks 2020)','lh','NDRTT.mat')}; % sourspot example 1 (tract))
sweetspot_weights=[1 ...
    1 ...
    -1
    -1];
useVTAvsEfields = 'vtas'; % set to 'efields' to use efields instead (both are always generated in step 1



%% run code
% Step 1: Calculate Pseudo VTAs
ea_mkdir('data');
if ~exist(['data',filesep,'generated_stimvols.mat'],'file') || reforce
    listout=ea_gen_pseudovtas('data',roidef,... % mask defining area covered by VTAs
        numvtas,... % number of vtas
        amp,...
        vta_jit); % VTA spread / jitter
    save(['data',filesep,'generated_stimvols.mat'],'listout');
else
    load(['data',filesep,'generated_stimvols.mat']);
end

% Step 2: Get Pseudo Improvements
if ~exist(['data',filesep,'improvements.mat'],'file') || reforce
    I=ea_gen_pseudoimprovements(listout.(useVTAvsEfields),...
        sweetspots,sweetspot_weights);
    
    I(isnan(I))=0; % neutral null
    I=ea_normal(I,1,1,1,'TRUE'); % gaussianize improvements
    save(['data',filesep,'improvements.mat'],'I');
else
    load(['data',filesep,'improvements.mat']);
end

% Export Groundtruth Sweetspot
if ~exist(fullfile('data',['groundtruth_sweetspot.nii']),'file') || reforce
    ea_gen_groudtruthsweetspot('data',sweetspots,sweetspot_weights);
end


% Step 3: Setup sweetspot / fiberfiltering explorers
improvements.groundtruth=I;
improvements.variations=repmat(I,1,1000);
improvements.variations=improvements.variations+...
    (randn(size(improvements.variations,1),size(improvements.variations,2))).*0.5;
%R=corr(improvements.variations);
%disp(['Average R = ',num2str(mean(R(:)))]);

resultfig = ea_mnifigure; % Create empty 3D viewer figure

% set up pseudo M struct:
M.pseudoM = 1; % Declare this is a pseudo-M struct, i.e. not a real lead group file
M.ROI.list = listout.efields'; % always use efields here (at least for sweetspot mapping
M.ROI.group=ones(length(M.ROI.list),1);
M.clinical.labels={'Groundtruth_Improvement'};
M.clinical.vars{1}=improvements.groundtruth;
for noisyI=1:size(improvements.variations,2)
    M.clinical.labels{noisyI+1}=['Noisy_Improvement_',num2str(noisyI)];
    M.clinical.vars{noisyI+1}=improvements.variations(:,noisyI);
end

M.guid='Pseudo_Sweetspot_Analysis'; % give your analysis a name
save(['data',filesep,'Analysis_Input_Data.mat'],'M'); % store data of analysis to file

ea_sweetspotexplorer(['data',filesep,'Analysis_Input_Data.mat'],resultfig);

% for Min Jae: instead of the above line run:
%ea_networkmappingexplorer(['data',filesep,'Analysis_Input_Data.mat'],resultfig);



