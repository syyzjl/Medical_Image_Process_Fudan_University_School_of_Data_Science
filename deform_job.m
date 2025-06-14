%-----------------------------------------------------------------------
% Job saved on 13-Jun-2025 22:33:27 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM25 (25.01.02)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

phi2 = 'C:\MIP_pj\data_code\data&code1.1\data\101799\y_101799.nii';
moving_path = 'C:\MIP_pj\data_code\data&code1.1\data\101799\101799.nii,1';
deformjob(phi2, moving_path);

function deformjob(phi, source)    
    % 重置图形系统
    set(0, 'DefaultFigureVisible', 'off');
    spm_figure('Clear', 'Graphics');
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');

    matlabbatch{1}.spm.spatial.normalise.write.subj.def = phi;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample =source;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                            78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

    nrun = 1; % 运行次数
    jobfile = {}; % 不使用外部job文件，直接使用内存中的matlabbatch
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(0, nrun);

    spm_jobman('run', matlabbatch);

    clear matlabbatch
    % 确保清除所有图形
    close all hidden;

end