function reg()
addpath('utils');
dirs = dir('.\data');
dirs(1:2) = [];
subjects = {dirs.name};
for i=1:length(subjects)
    %读取nii数据
    subject = subjects{i};
    subject_path = ['data\',subject];
    date_dir = dir(subject_path);
    mri_path = date_dir(~cell2mat({date_dir.isdir})).name;
    mri_path = [subject_path,'\', mri_path];
    date_dir = date_dir(cell2mat({date_dir.isdir}));
    date_dir(1:2) = [];
    
    for j=1:length(date_dir)
        moving_path = getAllFiles([subject_path,'\',date_dir(j).name]);                   
        moving_path = moving_path{1};
        % movingVolume = niftiread(moving_path);
        % moving_info = niftiinfo(moving_path);
        % movingVolume(isnan(movingVolume)) = 0;
        % fixedVolume = niftiread(mri_path);
        % fixed_info = niftiinfo(mri_path);
    
        %% 完成的配准代码_1的结果保存在这里
       % [movingRegisteredVolume, phi1]  = datscan_to_mri(movingVolume, ...
           % fixedVolume, moving_info,fixed_info);
        %将得到的配准结果保存为nii文件
        % temp = split(moving_path, '\');
        % temp_path = ['.\data\',temp{2},'\',temp{3},'\r',temp{4}];

        % niftiwrite(movingRegisteredVolume, temp_path, moving_info);
        coregister_job(moving_path, mri_path);

    end

    date_dir = dir(subject_path);
    mri_path = date_dir(~cell2mat({date_dir.isdir})).name;
    mri_path_ = [subject_path,'\', mri_path];
    date_dir = date_dir(cell2mat({date_dir.isdir}));
    date_dir(1:2) = [];
    
    tmp_path = '.\tpm\TPM.nii';
    % movingVolume = niftiread(mri_path);
    % moving_info = niftiinfo(mri_path);
    % movingVolume(isnan(movingVolume)) = 0;
    % fixedVolume = niftiread(tmp_path);
    % fixed_info = niftiinfo(tmp_path);
    %% 完成的配准代码_2的结果保存在这里
    normalize_(mri_path_, tmp_path);
     
    phi2  = [subject_path,'\y_', mri_path];
    

    date_dir = dir(subject_path);
    date_dir = date_dir(cell2mat({date_dir.isdir}));
    date_dir(1:2) = [];

    %% 对每个datscan作用这个变换
    for j=1:length(date_dir)
        moving_path = getAllFiles([subject_path,'\',date_dir(j).name]);
        moving_path = moving_path{1};
        % movingVolume = niftiread(moving_path);
        % moving_info = niftiinfo(moving_path);
        % movingVolume(isnan(movingVolume)) = 0;
        %% 对每个datscan作用这个变换
        % movingRegisteredVolume = transform2(movingVolume, phi1, phi2);
        %将得到的配准结果保存为rr为前缀的nii文件
        % temp = split(moving_path, '\');
        % temp_path = ['.\data\',temp{2},'\',temp{3},'\r',temp{4}];
        %% 可以保持PixelDimensions和ImageSize与moving_info中的值相同,或者根据配准代码_2修改moving_info对应的值
        % niftiwrite(movingRegisteredVolume, temp_path, moving_info);
        deform_job(phi2, moving_path);
    end



    
end



% 强度归一化

dirs = dir('.\data');
dirs(1:2) = [];
dirs = {dirs.name};

mask = niftiread('.\mask\roccipital.nii');
mask = double(mask);
mask_num = sum(mask,'all');

for i=1:length(dirs)
    dates = dir(['.\data\',dirs{i}]);
    dates = dates(cell2mat({dates.isdir}));
    dates(1:2) = [];
    dates = {dates.name};
    for k=1:length(dates)
        all_files = getAllFiles(['.\data\',dirs{i},'\',dates{k}]);
        image = all_files(cellfun(@(x)contains(x, '\w'), all_files));
        image = image{1};
        [filepath,name,ext] = fileparts(image);
        output_path = [filepath,'\n',name,ext];
        V = spm_vol([image,',1']);
        Y = spm_read_vols(V);
        temp = Y .* mask;
        const = nansum(temp,'all') / mask_num; %#ok<*NANSUM> 
        Y(~isnan(Y)) = (Y(~isnan(Y))-const) / const;
        V.dt = [16, 0];
        V.fname = output_path;
        spm_write_vol(V,Y);      
    end
end




%% 生成效果参考图片_version2
addpath('utils');


dirs = dir('.\data');
dirs(1:2) = [];
subjects = {dirs.name};
inds = 24:47;

caudate_margin = niftiread('.\mask\Margin Caudate.nii');
putamen_margin = niftiread('.\mask\Margin Putamen.nii');


for i=1:length(subjects)
    all_files = getAllFiles(['.\data\',subjects{i}]);
    mri_idx = cellfun(@(x)contains(x, '\w')&(~contains(x, 'DaTSCAN')), all_files);
    dat_idx = cellfun(@(x)contains(x, '\wr'), all_files);
    mri_path = all_files{mri_idx};
    V_mri = double(niftiread(mri_path));
    V_mri(isnan(V_mri)) = 0;
    V_mri = (V_mri - min(V_mri, [], 'all')) / (max(V_mri, [], 'all') - min(V_mri, [], 'all'));
    dat_path = all_files(dat_idx);
    for j=1:length(dat_path)
        figure(1);
        set(gcf, 'unit', 'centimeters', 'position', [3 3 30 11]);
        V_dat = double(niftiread(dat_path{j}));
        V_dat(isnan(V_dat)) = 0;
        V_dat = (V_dat - min(V_dat, [], 'all')) / (max(V_dat, [], 'all') - min(V_dat, [], 'all'));
        I_checkboard = zeros(89,105,length(inds));
        for k=1:length(inds)
            temp1 = padarray(V_mri(:,:,inds(k)), [5,5]);
            temp2 = padarray(V_dat(:,:,inds(k)), [5,5]);
            I_checkboard(:,:,k) = cbimage(temp1, temp2, 4);
        end
        subplot(1,3,1);
        montage(I_checkboard, 'size', [6, 4]);
        
        I_checkboard = zeros(89,105,3,length(inds));
        for k=1:length(inds)
            temp = repmat(padarray(V_dat(:,:,inds(k)), [5,5]),[1,1,3]);
            caudate_mask = repmat(padarray(caudate_margin(:,:,inds(k)), [5,5]),[1,1,3]);
            caudate_mask(:,:,3) = 0;
            putamen_mask = repmat(padarray(putamen_margin(:,:,inds(k)), [5,5]),[1,1,3]);
            putamen_mask(:,:,1) = 0;
            mask = imadd(caudate_mask, putamen_mask);
            I_checkboard(:, :, :, k) = imadd(double(mask)*0.3, temp);
        end
        subplot(1,3,2);
        montage(I_checkboard,'size',[6,4]);
        
        I_checkboard = zeros(89,105,length(inds));
        for k=1:length(inds)
            I_checkboard(:, :, k) = padarray(V_dat(:,:,inds(k)), [5,5]);
        end
        subplot(1,3,3);
        montage(I_checkboard,'size',[6,4]);

        temp = split(dat_path{j}, '\');
        sgtitle([temp{3},': ',temp{4}]);
        saveas(gcf,['./check_results_v2/',temp{3},'-',temp{4},'.jpg']);
    end
end

end


% function [movingRegisteredVolume, phi] = datscan_to_mri(movingVolume, fixedVolume, moving_info, fixed_info)
    % Rigid registration of DaTscan to MRI using SPM12
    % Input:
    %   movingVolume - DaTscan volume to register (3D matrix)
    %   fixedVolume  - Reference MRI volume (3D matrix)
    %   moving_info  - Moving volume metadata
    %   fixed_info   - Fixed volume metadata
    % Output:
    %   movingRegisteredVolume - Registered DaTscan volume
    %   phi - Transformation parameters
    
    
% end

function coregister_job(source, template)

    % 重置图形系统
    set(0, 'DefaultFigureVisible', 'off');
    spm_figure('Clear', 'Graphics');
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');

    source_file = source;      % 待配准的图像
    template_file = template;  % 模板图像

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
    clear matlabbatch
    % 确保清除所有图形
    close all hidden;
end



%% 需要完成的配准代码_2_获得mri到tpm的矩阵
%function phi = mri_to_tpm(movingVolume, fixedVolume, moving_info, fixed_info)
    % Non-rigid registration of MRI to template using SPM12
    % Input:
    %   movingVolume - MRI volume to register (3D matrix)
    %   fixedVolume  - Template volume (3D matrix)
    %   moving_info  - Moving volume metadata
    %   fixed_info   - Fixed volume metadata
    % Output:
    %   phi - Transformation parameters
    
    
% end

function normalize_(source, template)
    % 初始化SPM
    % 重置图形系统
    set(0, 'DefaultFigureVisible', 'off');
    spm_figure('Clear', 'Graphics');
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');
    
    % 定义文件路径
    source_file = source;  % 待配准的图像
    template_file = template;        % 模板图像
    
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

    clear matlabbatch
    % 确保清除所有图形
    close all hidden;
end


% function movingVolume1 = transform2(movingVolume, phi1, phi2)
    % Apply composite transformation to volume
    % Input:
    %   movingVolume - Volume to transform
    %   phi1 - First transformation (rigid)
    %   phi2 - Second transformation (non-rigid)
    % Output:
    %   movingVolume1 - Transformed volume
    
    % Create temporary file
    
%end

function deform_job(phi, source)    
    % 重置图形系统
    set(0, 'DefaultFigureVisible', 'off');
    spm_figure('Clear', 'Graphics');
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');

    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {phi};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample ={source};
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