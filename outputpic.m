function outputpic()
    addpath('utils');


    dirs = dir('.\data');
    dirs(1:2) = [];
    subjects = {dirs.name};
    inds = 24:47;

    caudate_margin = niftiread('.\mask\Margin Caudate.nii');
    putamen_margin = niftiread('.\mask\Margin Putamen.nii');


    for i=1:length(subjects)
        all_files = getAllFiles(['.\data\',subjects{i}]);
        % mri_idx = cellfun(@(x)contains(x, '\w')&(~contains(x, 'DaTSCAN')), all_files);
        mri_idx = cellfun(@(x) ~isempty(regexp(x, 'w\d{4,}', 'once')), all_files);
        % dat_idx = cellfun(@(x)contains(x, '\nw'), all_files);
        dat_idx = cellfun(@(x) contains(x, 'nw'), all_files);
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