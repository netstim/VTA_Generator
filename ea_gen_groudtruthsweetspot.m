function ea_gen_groudtruthsweetspot(outptFolder,sweetspots,weights)

if ~exist('weights','var')
    weights=ones(length(sweetspots),1);
end

% gather sweetspots:
for sws=1:length(sweetspots)
    [pth,fn,ext]=fileparts(sweetspots{sws});
    switch ext
        case {'.nii','.gz'} % roi / nucleus
            roi{sws}=ea_load_nii(sweetspots{sws});
        case {'.mat'} % fibertract
            ea_ftr2nii(sweetspots{sws}, [ea_space,'bb.nii'], fullfile(ea_getleadtempdir,['roi_',num2str(sws),'.nii']));
            roi{sws}=ea_load_nii(fullfile(ea_getleadtempdir,['roi_',num2str(sws),'.nii']));
    end
    roi{sws}.img=roi{sws}.img.*weights(sws);
    roi{sws}.fname=fullfile(ea_getleadtempdir,['roi_',num2str(sws),'.nii']);
    ea_write_nii(roi{sws});
    roilist{sws}=roi{sws}.fname;
end

% generate compound image:
matlabbatch{1}.spm.util.imcalc.input = [{[ea_space,'bb.nii']};
    roilist'];
matlabbatch{1}.spm.util.imcalc.output = fullfile(outptFolder,['groundtruth_sweetspot.nii']);
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = 'nansum(X(2:end,:))';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{1}.spm.util.imcalc.options.mask = -1;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',{matlabbatch});
ea_crop_nii(fullfile(outptFolder,['groundtruth_sweetspot.nii']));





