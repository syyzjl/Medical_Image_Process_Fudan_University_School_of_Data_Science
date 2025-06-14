function normalize_()
    % 初始化SPM
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');
    
    % 定义文件路径
    source_file = 'C:\MIP_pj\data_code\data&code1.1\data\40709\40709.nii';  % 待配准的图像
    template_file = 'C:\MIP_pj\data_code\data&code1.1\tpm\TPM.nii';        % 模板图像
    
    % 检查文件是否存在
    if ~exist(source_file, 'file')
        error('源文件不存在: %s', source_file);
    end
    if ~exist(template_file, 'file')
        error('模板文件不存在: %s', template_file);
    end
    
    % =============================================
    % 第一部分：设置批处理参数
    % =============================================
    matlabbatch{1}.spm.spatial.normalise.est.subj.vol = {source_file};
    matlabbatch{1}.spm.spatial.normalise.est.eoptions.biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.normalise.est.eoptions.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.normalise.est.eoptions.tpm = {template_file};
    matlabbatch{1}.spm.spatial.normalise.est.eoptions.affreg = 'mni';
    matlabbatch{1}.spm.spatial.normalise.est.eoptions.reg = [0 0 0.1 0.01 0.04];
    matlabbatch{1}.spm.spatial.normalise.est.eoptions.fwhm = 0;
    matlabbatch{1}.spm.spatial.normalise.est.eoptions.samp = 3;
    
    % =============================================
    % 第二部分：运行批处理
    % =============================================
    nrun = 1; % 运行次数
    jobfile = {}; % 不使用外部job文件，直接使用内存中的matlabbatch
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(0, nrun);
    
    % 运行批处理
    try
        spm_jobman('run', matlabbatch);
        disp('标准化处理完成！');
        
        % 检查输出文件
        [filepath, name, ext] = fileparts(source_file);
        output_file = fullfile(filepath, ['y_' name ext]);
        if exist(output_file, 'file')
            disp(['变形场已生成: ' output_file]);
        else
            warning('未找到预期的变形场输出文件');
        end
    catch ME
        disp('处理过程中出现错误:');
        disp(ME.message);
        rethrow(ME);
    end
end