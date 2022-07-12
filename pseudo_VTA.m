%% Pseudo-VTA Generator
% Min Jae Kim (mkim61@bwh.harvard.edu)
% Last Modified: 07/11/2022
%% INPUT (MUST IN IN THE SAME DIRECTOY PATH):
% (1) standard e-field template (Medtronic 3389)
% (2) electrode model (Medtronic 3389)
% (3) ea_stats.mat template 
% (4) stimparameters.mat template 
% (5) ROI Input (e.g STN)
%% OUTPUT: 
% (1) VTA.nii (Binarized VTA image)
% (2) VTA.mat 
% (3) VTA_Efield.nii (Efield image)
% (5) VTA_Efield_Gaussian (Gaussian-filtered Efield Image)
%% User Input Parameters 
stim_side="lh"; %choices: "lh" or "rh" (left or right STN volumes) 
ampi=2.5; %Amplitude (mA)  
jit_VTA_loc=0.5; %variance of VTA seed location from seed (MNI)
jit_VTA_amp=0.5; %variance of VTA seed radius by current amplitude (mA) 
num_vtas=1000; % numberof VTAs generated (MAX: 4000 VTAS)
%% Loading Necessary Files
addpath(genpath("/Users/mkim197/Desktop/Lead_NEW")); %Add path to Lead-DBS
load('/Users/mkim197/Desktop/Lead_NEW/lead/templates/standard_efields/standard_efield_medtronic_3389.mat'); %standard_efield template for Medtronic 3389
load('/Users/mkim197/Desktop/Lead_NEW/lead/templates/electrode_models/medtronic_3389.mat'); %electrode model for Medtronic 3389
load('ea_stats.mat') %template form of ea_stats.mat 
load('stimparameters.mat') %template form of stimparameters.mat 
nii=ea_load_nii(char("rre_"+stim_side+"_STN.nii")); %input ROI_name: "re_rh_STN.nii"
%% 
max_val=max(nii.img(:));
min_val=min(nii.img(:));
[X,Y,Z]=ind2sub(size(nii.img),find(nii.img(:)>min_val));
XYZvox=[X,Y,Z,ones(size(X,1),1)];
XYZmm=nii.mat*XYZvox';
XYZmm=XYZmm(1:3,:)';

%
if stim_side=="rh"
    side =1;
    S.Rs1.amp=ampi;
elseif stim_side=="lh"
    side =2;
    S.Ls1.amp=ampi;
end 

for ggg=1:num_vtas %num_vtas
%% User Input
idxi=randi(size(XYZmm,1));
coor=XYZmm(idxi,:); 
coor=[normrnd(coor(1),jit_VTA_loc),normrnd(coor(2),jit_VTA_loc),normrnd(coor(3),jit_VTA_loc)];
ampi=normrnd(ampi,jit_VTA_amp); %Amplitude (mA) + jitter
%
% 
options=ea_getptopts(pwd);
useSI = 1;
acoords=S.activecontacts; %varargin{1};
S=S;
options=options;
stimname='SAMPLE';

conductivity = options.prefs.machine.vatsettings.fastfield_cb;  % 0.1;
thresh = options.prefs.machine.vatsettings.fastfield_ethresh; % 0.2;

if useSI
    thresh=thresh.*(10^3);
end

if ~any(S.activecontacts{side}) % empty VAT, no active contacts.
    ofv.vertices=[0,0,0
        0,0,0
        0,0,0];
    ofv.faces=[1,2,3];
    varargout{1}=ofv;
    varargout{2}=0;
    return
end

elstruct=ea_stats.electrodes; 
options=options; 
options=ea_resolve_elspec(options);
elspec=options.elspec; 
options.usediffusion=0;
coords=S.activecontacts{side}; 
%setappdata(resultfig,'elstruct',elstruct);

switch side
    case 1
        sidec='R';
        cnts={'k0','k1','k2','k3','k4','k5','k6','k7'};
    case 2
        sidec='L';
        cnts={'k8','k9','k10','k11','k12','k13','k14','k15'};
end

if ~isfield(S, 'sources')
    S.sources=1:4;
end

Electrode_type = elspec.matfname;

Efield_all=zeros(100,100,100);
for source=S.sources
    stimsource=S.([sidec,'s',num2str(source)]);
    % constvol is 1 for constant voltage and 0 for constant current.
    amp1 = stimsource.amp;
    if amp1>0
        %load('/Users/mkim197/Desktop/Lead_NEW/lead/templates/standard_efields/standard_efield_medtronic_3389.mat');
        %load([ea_getearoot,'templates',filesep,'standard_efields' filesep 'standard_efield_' Electrode_type '.mat']);
        count1=1;
        for cnt=1:length(cnts)
            perc(cnt) = stimsource.(cnts{cnt}).perc; %percentage of stim distribution across contacts
            if perc(cnt)>0
                Im(count1)=stimsource.(cnts{cnt}).imp; %impedance for stim distribution across contacts
                count1=count1+1;
            end
        end

        constvol=stimsource.va==1;
        if constvol %check if it is voltage stim 
            amp_mode = 'V';
            impedence = mean(Im)*1000;
        else        %check if it is amp stim 
            amp_mode = 'mA';
            impedence = [];
        end

        [Efield2] = get_efield(perc,standard_efield,amp1,conductivity,amp_mode,impedence);
        Efield_all=Efield_all+Efield2; %gaussian addition of e-field generated across each contacts
    end
end

Efield = Efield_all;

electrode_patient = elstruct;
%% 

if side==1
    electrode_patient.markers(1).head=coor+[-2.1715 -4.24653 0.02761];
    electrode_patient.markers(1).tail= electrode_patient.markers(1).head+[3.23    3.40    9.4]; %[3.2268    3.9788    9.2692];
    electrode_patient.markers(1).x = electrode_patient.markers(1).head.*[1.0649    1.0044    1.0107]; 
    electrode_patient.markers(1).y = electrode_patient.markers(1).head.*[1.0000    0.9663    1.0151];

elseif side ==2
    electrode_patient.markers(2).head=coor+[2.1806 -1.1825 -4.7097];  
    electrode_patient.markers(2).tail=electrode_patient.markers(2).head+[-3.27    4.3658    8];  
    electrode_patient.markers(2).x = electrode_patient.markers(2).head.*[0.9354    0.9949    0.9895];
    electrode_patient.markers(2).y = electrode_patient.markers(2).head.*[1.0000    0.9654    1.0163];
end
%% 
[trans_mat,~,xg,yg,zg] = get_trans_mat(electrode,electrode_patient,grid_vec,side);

gv=grid_vec;

ea_dispt('Creating nifti header for export...');
% create nifti
res=100;
chun1=randperm(res); chun2=randperm(res); chun3=randperm(res);
Vvat.mat=mldivide([(chun1);(chun2);(chun3);ones(1,res)]',[gv{1}(chun1);gv{2}(chun2);gv{3}(chun3);ones(1,res)]')';
Vvat.mat = trans_mat * Vvat.mat;
Vvat.dim=[res,res,res];
Vvat.dt=[4,0];
Vvat.n=[1 1];
Vvat.descrip='lead dbs - vat';
if ~exist([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options)],'file')
    mkdir([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options)]);
end

eeg = Efield;
eg=eeg;
eg=eg>thresh;
% binary e-field - "vat"

neeg=eeg;
neeg(~eg)=nan;

neeg(neeg>0)=ea_normal(neeg(neeg>0),1,0,' ',0,1,'TRUE');%
% normalized e-field (zscored).
neeg(~isnan(neeg))=neeg(~isnan(neeg))-min(neeg(~isnan(neeg)));
neeg(~isnan(neeg))=neeg(~isnan(neeg))/sum(neeg(~isnan(neeg))); % 0-1 distributed.

% [xg,yg,zg] = meshgrid(gv{1},gv{2},gv{3});

% eg=smooth3(eg,'gaussian',[25 25 25]);
ea_dispt('Calculating volume...');

for t=1:3
    spacing(t) = grid_vec{t}(2)-grid_vec{t}(1);
end
vatvolume=sum(eg(:))*spacing(1)*spacing(2)*spacing(3); % returns volume of vat in mm^3
S.volume(side)=vatvolume;

ea_dispt('Writing files...');

% determine stimulation name:
if ~exist([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),filesep,stimname],'file')
    mkdir([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),filesep,stimname]);
end

switch side
    case 1
        vat_tit=char("vat_right_"+string(ggg)+".nii");
        vat_efield=char("vat_efield_right_"+string(ggg)+".nii");
        vat_gaus=char("vat_efield_gauss_right_"+string(ggg)+".nii");
        Vvat.fname=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),filesep,stimname,filesep,vat_tit];
        Vvate=Vvat; Vvatne=Vvat;
        Vvate.fname=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),filesep,stimname,filesep,vat_efield];
        Vvatne.fname=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),filesep,stimname,filesep,vat_gaus];
    case 2
        vat_tit=char("vat_left_"+string(ggg)+".nii");
        vat_efield=char("vat_efield_left_"+string(ggg)+".nii");
        vat_gaus=char("vat_efield_gauss_left_"+string(ggg)+".nii");
        Vvat.fname=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),filesep,stimname,filesep,vat_tit];
        Vvate=Vvat; Vvatne=Vvat;
        Vvate.fname=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),filesep,stimname,filesep,vat_efield];
        Vvatne.fname=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),filesep,stimname,filesep,vat_gaus];
end

% save(stimfile,'S');
ea_savestimulation(S,options);
% setappdata(lgfigure,'curS',S);

%spm_write_vol(Vvat,flipdim(eg,3));

Vvate.img=eeg; %permute(eeg,[2,1,3]);
Vvate.dt=[16,0];
ea_write_nii(Vvate);

Vvatne.img=neeg; %permute(neeg,[2,1,3]);
ea_write_nii(Vvatne);

Vvat.img=eg; %permute(eg,[1,2,3]);
ea_write_nii(Vvat);

ea_dispt('Calculating isosurface to display...');
% vatfv=isosurface(xg,yg,zg,permute(Vvat.img,[2,1,3]),0.75);
vatfv=isosurface(xg,yg,zg,Vvat.img,0.75);

% caps=isocaps(xg,yg,zg,permute(Vvat.img,[2,1,3]),0.5);
caps=isocaps(xg,yg,zg,Vvat.img,0.5);

vatfv.faces=[vatfv.faces;caps.faces+size(vatfv.vertices,1)];
vatfv.vertices=[vatfv.vertices;caps.vertices];

try
    vatfv=ea_smoothpatch(vatfv,1,35);
catch
    try
        cd([ea_getearoot,'ext_libs',filesep,'smoothpatch']);
        mex ea_smoothpatch_curvature_double.c -v
        mex ea_smoothpatch_inversedistance_double.c -v
        mex ea_vertex_neighbours_double.c -v
        vatfv=ea_smoothpatch(vatfv);
    catch
        warndlg('Patch could not be smoothed. Please supply a compatible Matlab compiler to smooth VTAs.');
    end
end

% new save by Till to save VAT and quiver in seperate .mat-file for quick
% visualization
switch side
    case 1
        vat_fv_tit=char("vat_right_"+string(ggg)+".mat");
        vatfvname=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),filesep,stimname,filesep,vat_fv_tit];
    case 2
        vat_fv_tit=char("vat_left_"+string(ggg)+".mat");
        vatfvname=[options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),filesep,stimname,filesep,vat_fv_tit];
end

save(vatfvname,'vatfv','vatvolume');

%% new vta.nii save, filled and eroded/dilated by 3 voxels.
Vvat.img=imfill(Vvat.img,'holes');
SE = strel('sphere',3);
Vvat.img = imerode(Vvat.img,SE);
Vvat.img = imdilate(Vvat.img,SE);
ea_write_nii(Vvat);

varargout{1}=vatfv;
varargout{2}=vatvolume;
ea_dispt('');

clc

end 
