%-----------------------------------------------------------------------
% Job saved on 13-Jun-2025 21:29:10 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM25 (25.01.02)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
function coregister_job()
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');

    source_file = 'C:\MIP_pj\data_code\data&code1.1\data\101799\1\1.nii';      % 待配准的图像
    template_file = 'C:\MIP_pj\data_code\data&code1.1\data\101799\101799.nii';  % 模板图像

    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {template_file};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {source_file};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    nrun = 1; % 运行次数
    jobfile = {}; % 不使用外部job文件，直接使用内存中的matlabbatch
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(0, nrun);
    
    spm_jobman('run', matlabbatch);

end