function listout=ea_gen_pseudovtas(varargin)
% Code to generate a set of fake VTAs
% usage:
% listout=ea_gen_pseudovtas(outptFolder,maskFn,numvtas,amp,elmodel)
%   outputfolder = where to store the outputs
%   mask = nifti defining centers (e.g. STN)
%   numvta = how many vtas
%   amp.mean = mean amplitude
%   amp.jit = amplitude jitter
%   OPTIONAL args:
%   elmodel = medtronic_3389 (default)
%   jitVTA = 0.5 (default)
%   OUTPUT arg:
%   listout.vtas / listout.efields = cell strings of files written out
%   EXAMPLE:
%      
% amp.mean=2.5;
% amp.jit=2;
% ea_gen_pseudovtas(fullfile(ea_space([],'atlases'),'DISTAL Nano (Ewert 2017)','lh','STN.nii.gz'),...
%     1000,{fullfile(ea_space([],'atlases'),'TOR-signPD (Boutet 2021)','lh','tremor_hotspot.nii.gz'),...
%     fullfile(ea_space([],'atlases'),'DBS Tractography Atlas (Middlebrooks 2020)','lh','DRTT.mat')});
%     
%
%

if nargin<5
    ea_error('Not enough input arguments provided.');
end

options=ea_getptopts(pwd);
outptFolder=varargin{1};
maskFn=varargin{2};
numvtas=varargin{3};
amp=varargin{4};
if nargin>4
    jitVTA=varargin{5};
else
    jitVTA=0.5;
end

if nargin>5
    elmodel=varargin{6};
else
    elmodel='medtronic_3389';
end




% define seed centroids based on mask:
mask=ea_load_nii(maskFn); 
[xx,yy,zz]=ind2sub(size(mask.img),find(mask.img(:)>0));
XYZ=[xx,yy,zz,ones(size(xx,1),1)]';
centroidMm=mask.mat*XYZ;
centroidMm=centroidMm(1:3,:)';
% figure, plot3(XYZmm(:,1),XYZmm(:,2),XYZmm(:,3),'r.')
numCentroids=size(centroidMm,1);
if numCentroids>numvtas
    centroidMm=centroidMm(round(1:(numCentroids/numvtas):numCentroids),:);
    centroidMm=centroidMm(1:numvtas,:); % rounding errors..
    numCentroids=size(centroidMm,1);
end

vtasPerCentroid=upper(numvtas/numCentroids);

% load fastfield standard efield:
ef=load(fullfile(ea_getearoot,'templates','standard_efields',['standard_efield_',elmodel,'.mat']));
mdl=load(fullfile(ea_getearoot,'templates','electrode_models',[elmodel,'.mat']));

% define VTA threshold:
thresh = options.prefs.machine.vatsettings.fastfield_ethresh.*(10^3); % useSI



% iterate through all centroids to generate #vtasPerCentroid in each spot:
cnt=1;
ea_dispercent(0,'Writing out efields & vtas')
for cent=1:numCentroids
    centroid=centroidMm(cent,:); % current centroid
    for vta=1:vtasPerCentroid
        coord=[normrnd(centroid(1),jitVTA),normrnd(centroid(2),jitVTA),normrnd(centroid(3),jitVTA)];
        thisamp=-1; % make sure no negative amplitude gets created by jitter
        while thisamp<0
            thisamp=normrnd(amp.mean,amp.jit);
        end
        
        [Efield] = ea_get_efield([100;0;0;0],ef.standard_efield,thisamp,...
            options.prefs.machine.vatsettings.fastfield_cb,...
            'mA',[]);

        res=size(Efield,1); % assuming isotropic resolution of the fastfield output
        chun1=randperm(res); chun2=randperm(res); chun3=randperm(res);
        nii = ea_getdefnii;
        nii.mat=mldivide([(chun1);(chun2);(chun3);ones(1,res)]',[ef.grid_vec{1}(chun1);ef.grid_vec{2}(chun2);ef.grid_vec{3}(chun3);ones(1,res)]')';

        nii.mat(3,4)=nii.mat(3,4)-2; % offset t o make center of efield 0,0,0
        coord=nii.mat(1:3,4)'+coord;
        nii.mat(1:3,4)=coord';
        nii.img=Efield;
        nii.fname=fullfile(outptFolder,['efield_',num2str(cnt),'.nii']);
        listout.efields{cnt}=nii.fname;
        ea_write_nii(nii);
        nii.img=nii.img>thresh;
        nii.fname=fullfile(outptFolder,['vta_',num2str(cnt),'.nii']);
        ea_write_nii(nii);
        listout.vtas{cnt}=nii.fname;
        cnt=cnt+1;
    end
    ea_dispercent(cent/numCentroids);
end
ea_dispercent(1,'end');


% generate N-image
matlabbatch{1}.spm.util.imcalc.input = [{[ea_space,'bb.nii']};
    listout.vtas'];
matlabbatch{1}.spm.util.imcalc.output = fullfile(outptFolder,['nimage.nii']);
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = 'nansum(X(2:end,:))';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{1}.spm.util.imcalc.options.mask = -1;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',{matlabbatch});
ea_crop_nii(fullfile(outptFolder,['nimage.nii']));

matlabbatch{1}.spm.util.imcalc.input = [{[ea_space,'bb.nii']};
    listout.efields'];
matlabbatch{1}.spm.util.imcalc.output = fullfile(outptFolder,['nimage_efield.nii']);
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = 'nansum(X(2:end,:))';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{1}.spm.util.imcalc.options.mask = -1;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',{matlabbatch});
ea_crop_nii(fullfile(outptFolder,['nimage_efield.nii']));

save(fullfile(outptFolder,['generated_stimvols.mat']),'listout');





function nii=ea_getdefnii
        nii.fname = '';
        nii.dim = [100 100 100];
        nii.dt = [4 0];
        nii.pinfo = [1;             0;             352];
        nii.mat = [1 0 0 -50
            0 1 0 -50
            0 0 1 -46
            0 0 0 1];
        nii.n = [1 1];
        nii.descrip = '';
        nii.voxsize = [1 1 1];
        nii.volnum = 1;




