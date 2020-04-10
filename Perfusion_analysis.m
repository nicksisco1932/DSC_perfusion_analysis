%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The whole idea behind this script is to make it easy to load the MRI data
% in order to analyze it.

% This script asks the user to first identify a base directory where you
% all of the patient directories are located. Then it asks which folders
% you want to process.

% You have the option of clicking one directory or selecting all of them,
% or even CTRL+ selecting directories asynchronously. It is my hope that
% this will prevent the need from having many different scripts to analyze
% the data and make it so there is less errors in the process.

% Created by Nicholas J. Sisco, Ph.D. please send questions to
% nicholas dot sisco at barrowneuro dot org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;close all;clear
disp('Select the base folder')
temp = uigetdir;
path = dir(temp);

cd(path(1).folder)
base_dir = pwd;
clc
disp('Select all of the directories you wish to process')
lots_of_dirs = uigetfile_n_dir();

%
lots_of_dirs{1};
tmp = table(lots_of_dirs);
tmp.lots_of_dirs{:,1};

N = size(lots_of_dirs);
%
for i = 1:N(2)
% for i = 3
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    cd(tmp.lots_of_dirs{:,i});clear T2w_FLAIR T1w_PREC T1w_POSTC DSC
    a = strsplit(tmp.lots_of_dirs{:,i},'/');disp(['Processing ' a{end}]);
    disp(['Processing ' num2str(i) ' of ' num2str(N(2)) ' selected directories'])
    
    disp('Data Load')
    save_file_flag = 0; % easier than typing it for every function.
    try
        [T2w_FLAIR.Data, T2w_FLAIR.Spc] = t2FLAIR(i,tmp,save_file_flag );
        tempI = make_nii(imrotate3(T2w_FLAIR.Data,-90,[0 0 1]),T2w_FLAIR.Spc,[0 0 0]);
        save_nii(tempI, [pwd  '/T2w_FLAIR.nii'], 0);
    catch
        warning(['Check on ' a{end} 'T2w_FLAIR']);
    end
    try
        [T1w_pre.Data, T1w_pre.Spc] = t1PRE(i,tmp,save_file_flag );
        tempI = make_nii(imrotate3(T1w_pre.Data,-90,[0 0 1]),T1w_pre.Spc,[0 0 0]);
        save_nii(tempI, [pwd  '/T1w_pre.nii'], 0);
    catch
        warning(['Check on ' a{end} ' T1w_PREC']);
    end
    try
        [T1w_post.Data, T1w_post.Spc] = t1POST(i,tmp,save_file_flag );
        tempI = make_nii(imrotate3(T1w_post.Data,-90,[0 0 1]),T1w_post.Spc,[0 0 0]);
        save_nii(tempI, [pwd  '/T1w_post.nii'], 0);
    catch
        warning(['Check on ' a{end} 'T1w_POSTC']);
    end
    try
        [DSC.Data, DSC.Spc] = SAGE(i,tmp,save_file_flag );
        tempI = make_nii(imrotate3(DSC.Data,-90,[0 0 1]),DSC.Spc,[0 0 0]);
        save_nii(tempI, [pwd  '/SAGE.nii'], 0);
    catch
        warning(['Check on ' a{end} 'DSC']);
    end

    
%     % Select directories from this list first. This is just something weird
%     % with the data.
%     T = table(['PT1319001';'PT1319003';'PT1319004';'PT1319005';'PT1319007';'PT1319009';...
%     'PT1319010';'PT1319011';'PT1319012';'PT1319015';'PT1319016';'PT1319017';'PT1319019';...
%     'PT1319030';'PT1319031';'PT1319032';'PT1319033';'PT1319034';'PT1319035';'PT1319036';...
%     'PT1319037';'PT1319038';'PT1319040';'PT1319051';'PT1319047';'PT1319048';'PT1319052';...
%     'PT1319053';'PT1319055';'PT1319059';'PT1319062';'PT1319063';'PT1319065';'PT1319068';...
%     'PT1319073';]);
%     T2 = table(['PT1319005';'PT1319006';'PT1319007';'PT1319023';'PT1319026'...
%     ;'PT1319027';'PT1319028';'PT1319030';'PT1319031';'PT1319032'...
%     ;'PT1319033';'PT1319034';'PT1319035';'PT1319036';'PT1319037'...
%     ;'PT1319038';'PT1319048';'PT1319051';'PT1319046']);
%     n = size(T);
%     
%     T3 = table(['PT1319003';'PT1319052';'PT1319053';'PT1319055';'PT1319059';...
%         'PT1319060';'PT1319062';'PT1319063';'PT1319064';'PT1319065';...
%         'PT1319066';'PT1319068';'PT1319073';'PT1319047';'PT1319048']); %for i = 1:14 [96 96 150 5 15]
% 
%     T4 = 'PT1319003';
%
% level = 1;
%         for j = [6 ]
%             if strcmp(tmp.lots_of_dirs{:,i},join([base_dir, '/' T.Var1(j,:)])) == 1
%                 disp('YAAS')
%                 try
%                     [DSC.Data,DSC.Parms,DSC.AIF,DSC.AIFMask,DSC.dR2s,DSC.Perfusion] = conv_CBV_CBF1(DSC.Data);close all
%                     tempI = make_nii(imrotate3(DSC.Data,-90,[0 0 1]),DSC.Spc,[0 0 0]);
%                     save_nii(tempI, [pwd  '/' a{end} '_SAGE_.nii'], 0);
%                 catch
%                     warning(['Check on ' T.Var1(j,:) ' level ' num2str(level)])
%                 
%                 end
%             end
%         end
% level = 2;
%         for j = [1 7 8 9 10 11 12 13 14 ]
%             if strcmp(tmp.lots_of_dirs{:,i},join([base_dir, '/' T.Var1(j,:)])) == 1
%                 disp('YAAS')
%                 try
%                     [DSC.Data,DSC.Parms,DSC.AIF,DSC.AIFMask,DSC.dR2s,DSC.Perfusion] = conv_CBV_CBF2(DSC.Data);close all
%                     tempI = make_nii(imrotate3(DSC.Data,-90,[0 0 1]),DSC.Spc,[0 0 0]);
%                     save_nii(tempI, [pwd  '/' a{end} '_SAGE_.nii'], 0);
%                 catch
%                     warning(['Check on ' T.Var1(j,:) ' level ' num2str(level)])
%                 
%                 end
%             end
%         end
% level = 3;
%         for j = 1:19
%             if strcmp(tmp.lots_of_dirs{:,i},join([base_dir, '/' T2.Var1(j,:)])) == 1
%                 disp('YAAS')
%                 try
%                     [DSC.Data,DSC.Parms,DSC.AIF,DSC.AIFMask,DSC.dR2s,DSC.Perfusion] = conv_CBV_CBF3(DSC.Data);
%                     close all
%                     tempI = make_nii(imrotate3(DSC.Data,-90,[0 0 1]),DSC.Spc,[0 0 0]);
%                     save_nii(tempI, [pwd  '/' a{end} '_SAGE_.nii'], 0);
%                 catch
%                     warning(['Check on ' T2.Var1(j,:) ' level ' num2str(level)])
%                 
%                 end
%             end
%         end
% level = 4;
%         for j = 1:15
%             if strcmp(tmp.lots_of_dirs{:,i},join([base_dir, '/' T3.Var1(j,:)])) == 1
%                 disp('YAAS')
%                 try
%                     [DSC.Data,DSC.Parms,DSC.AIF,DSC.AIFMask,DSC.dR2s,DSC.Perfusion] = conv_CBV_CBF4(DSC.Data);
%                     close all
%                     tempI = make_nii(imrotate3(DSC.Data,-90,[0 0 1]),DSC.Spc,[0 0 0]);
%                     save_nii(tempI, [pwd  '/' a{end} '_SAGE_.nii'], 0);
%                 catch
%                     warning(['Check on ' T3.Var1(j,:) ' level ' num2str(level)])
%                 
%                 end
%             end
%         end
%         
%     disp(['Data Load Completed for ' a{end}]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The registration is done better in Python and I have already written this
% script.
%     disp(['Register images for ' a{end}]);
%     
%     fsl = getenv('FSLDIR'); %set the path for FSL
%     unix([fsl '/bin/bet ' a{end} '_T1w_pre_.nii brainT1w_pre.nii -f 0.5 -g 0.4 -m'])
%     unix([fsl '/bin/bet ' a{end} '_T1w_post_.nii brainT1w_post.nii -f 0.5 -g 0.4 -m'])
%     unix([fsl '/bin/bet ' a{end} '_T2w_FLAIR_.nii brainT2w_FLAIR.nii -f 0.4 -g 0 -m'])
%     unix([fsl '/bin/bet ' a{end} '_SAGE_.nii brainSAGE.nii -f 0.4 -g 0 -m'])
%     
%     in = 'brainT1w_post.nii';
%     out = 'AffineReg12p_brainT1wpost2pre.nii';
%     ref = 'brainT1w_pre.nii';
%     omat = 'Reg_refT1wpre_inputT1wpost.txt';
%     opts = ' -bins 256 -cost normmi -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear ';
%     opts_files = [' -in ' in ' -ref ' ref ' -out ' out ' -omat ' omat];
%     unix([fsl '/bin/flirt' opts_files opts]);
%     disp(['Phase 1 Registration Done on ' a{end}])
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,Spc] = t1PRE(N,directories,save_flag)
disp('T1 PRE')
% info = dicominfo('AX_FLAIR/IM-0001-0001-0001.dcm');
% info.ImagePositionPatient;

tmp = directories;
data.DATA = zeros(512,512,170);

temp = dir();
n = size(temp);
for i = 1:n(1)
    if strcmp(fullfile(temp(i).name),'T1w_TFE_preC') ==1
        file_number = i;
    elseif strcmp(fullfile(temp(i).name),'WIP AX 3D T1') ==1
        file_number = i;
    elseif strcmp(fullfile(temp(i).name),'WIP AX T1 3D_TFE (WAND) 1_2') ==1
        file_number = i;
    elseif strcmp(fullfile(temp(i).name),'WIP AX T1 3D WAND') ==1
        file_number = i;
    end
end

data_dir = fullfile(tmp.lots_of_dirs{:,N},fullfile(temp(file_number).name));
sub_dir = dir([data_dir,'/Series*']);

if ~isempty(sub_dir)
    subby = 1;
else
    subby = 0;
end

switch subby
    case 1
        if length(sub_dir) > 1
            sub_dir = sub_dir(end); % this just selects the last directory of the Series because the first directory was used as a test.
        end
        dir_list = dir(fullfile(data_dir,sub_dir.name));
        for i = 1:length(dir_list)
            ind(i,:) = ~startsWith(dir_list(i).name,'.') && ~endsWith(dir_list(i).name,'.dcm') && startsWith(dir_list(i).name, 'Image');
        end
        
        if sum(ind) == 0
            for i = 1:length(dir_list)
                ind(i,:) = ~startsWith(dir_list(i).name,'.') && ~endsWith(dir_list(i).name,'.dcm') && startsWith(dir_list(i).name, 'Mag');
            end
        end
        
        if sum(ind) == 0 % if it still didn't get the files. 
            for i = 1:length(dir_list)
                ind(i,:) = ~startsWith(dir_list(i).name,'.') && endsWith(dir_list(i).name,'.dcm') && startsWith(dir_list(i).name, 'I');
            end
        end
        
        count = 0;
        for i = 1:length(ind)
            if ind(i) == 1
                if count == 0
                    info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
                    try
                        Spc = [info.PixelSpacing; info.SliceThickness]';
                    catch
                        Spc = [0.5078 0.5078 1]; % reslice it after saving
                    end
                    count = count + 1;
                end
            end
        end
        count = 0;
        if sum(ind) > 1
            temp = dicomread(fullfile(data_dir,sub_dir.name,dir_list(5).name));
            [x,y] = size(temp);
            data.DATA = zeros(x,y,sum(ind));clear temp x y
        end
        if sum(ind) == 1
            for i = 1:length(ind)
                if ind(i) == 1
%                     info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
%                     try
%                         Spc = [info.PixelSpacing; info.SliceThickness]';
%                     catch
%                         Spc = [0.5078 0.5078 1]; % reslice it after saving
%                     end
                    temp = dicomread(fullfile(data_dir,sub_dir.name,dir_list(i).name));
                    data.DATA = squeeze(temp);
                    disp(['Reading image number ' num2str(count) ' of ' num2str(170)])
                end
            end
        else
            for i = 1:length(ind)
                if ind(i) == 1
%                     info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
%                     try
%                         Spc = [info.PixelSpacing; info.SliceThickness]';
%                     catch
%                         Spc = [0.5078 0.5078 1]; % reslice it after saving
%                     end
                    count = count+1;
                    temp = dicomread(fullfile(data_dir,sub_dir.name,dir_list(i).name));
                    data.DATA(:,:,count) = squeeze(temp);
                    disp(['Reading image number ' num2str(count) ' of ' num2str(170)])
                end
            end
        end
    case 0
        dir_list = dir(data_dir);
        for i = 1:length(dir_list)
            ind(i,:) = ~startsWith(dir_list(i).name,'.') && endsWith(dir_list(i).name,'.dcm');
        end
        
        count = 0;
        for i = 1:length(ind)
            if ind(i) == 1
                if count == 0
                    info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
                    try
                        Spc = [info.PixelSpacing; info.SliceThickness]';
                    catch
                        Spc = [0.5078 0.5078 1]; % reslice it after saving
                    end
                    count = count + 1;
                end
            end
        end
        count = 0;
        for i = 1:length(ind)
            if ind(i) == 1
                count = count+1;
%                 info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
%                     try
%                         Spc = [info.PixelSpacing; info.SliceThickness]';
%                     catch
%                         Spc = [0.5078 0.5078 1]; % reslice it after saving
%                     end
                temp = dicomread(fullfile(data_dir,dir_list(i).name));
                data.DATA(:,:,count) = squeeze(temp);
                disp(['Reading image number ' num2str(count) ' of ' num2str(170)])
            end
        end
        
end

data = data.DATA;
switch save_flag
    case 1
        save(fullfile(tmp.lots_of_dirs{:,1},[date '_T1w_PRE.mat']),'data','-mat')
end

end

function [data,Spc] = t2FLAIR(N,directories,save_flag)
disp('T2 Flair')
% info = dicominfo('AX_FLAIR/IM-0001-0001-0001.dcm');
% info.ImagePositionPatient;

tmp = directories;
data.DATA = zeros(512,512,30);

temp = dir();
n = size(temp);
for i = 1:n(1)
    if strcmp(fullfile(temp(i).name),'AX_FLAIR') ==1
        file_number = i;
    elseif strcmp(fullfile(temp(i).name),'WIP AX FLAIR MVXD') ==1
        file_number = i;
    elseif strcmp(fullfile(temp(i).name),'WIP AX T2 FLAIR MVXD') ==1
        file_number = i;
    elseif strcmp(fullfile(temp(i).name),'WIP AX T2 MVXD') ==1
        file_number = i;
    end
end

data_dir = fullfile(tmp.lots_of_dirs{:,N},fullfile(temp(file_number).name));
sub_dir = dir([data_dir,'/Series*']);

if ~isempty(sub_dir)
    subby = 1;
else
    subby = 0;
end

switch subby
    case 1
        %         warning('Code this for a subdiretory')
        if length(sub_dir) > 1
            sub_dir = sub_dir(end); % this just selects the last directory of the Series because the first directory was used as a test.
        end
        dir_list = dir(fullfile(sub_dir.folder,sub_dir.name));
        for i = 1:length(dir_list)
            % This will open the dicom image that does not have the dcm at
            % the end of it, which is how the data is originally when
            % unzipped. This should handle single files of this name.
            ind(i,:) = ~startsWith(dir_list(i).name,'.') && ~endsWith(dir_list(i).name,'.dcm') && startsWith(dir_list(i).name, 'Image');
        end
        
        if sum(ind) == 0
            for i = 1:length(dir_list)
                % This will open the dicom image that does not have the dcm at
                % the end of it, which is how the data is originally when
                % unzipped. This should handle single files of this name.
                ind(i,:) = ~startsWith(dir_list(i).name,'.') && ~endsWith(dir_list(i).name,'.dcm') && startsWith(dir_list(i).name, 'Mag');
            end
        end
        
        if sum(ind) == 0 % if it still didn't get the files. 
            for i = 1:length(dir_list)
                % This will open the dicom image that does not have the dcm at
                % the end of it, which is how the data is originally when
                % unzipped. This should handle single files of this name.
                ind(i,:) = ~startsWith(dir_list(i).name,'.') && endsWith(dir_list(i).name,'.dcm') && startsWith(dir_list(i).name, 'I');
            end
        end
        
        count = 0;
        for i = 1:length(ind)
            if ind(i) == 1
                if count == 0
                    info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
                    try
                        Spc = [info.PixelSpacing; info.SliceThickness]';
                    catch
                        Spc = [0.5078 0.5078 1]; % reslice it after saving
                    end
                    count = count + 1;
                end
            end
        end
        count = 0;
        if sum(ind) > 1
            temp = dicomread(fullfile(data_dir,sub_dir.name,dir_list(5).name));
            [x,y] = size(temp);
            data.DATA = zeros(x,y,sum(ind));clear temp x y
        end
        
        if sum(ind) == 1
            for i = 1:length(ind)
                if ind(i) == 1
%                     info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
%                     try
%                         Spc = [info.PixelSpacing; info.SliceThickness]';
%                     catch
%                         Spc = [0.5078 0.5078 1]; % reslice it after saving
%                     end
                    temp = dicomread(fullfile(data_dir,sub_dir.name,dir_list(i).name));
                    data.DATA = squeeze(temp);
                    disp(['Reading image number ' num2str(count) ' of ' num2str(30)])
                end
            end
        else
            for i = 1:length(ind)
                if ind(i) == 1
%                     info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
%                     try
%                         Spc = [info.PixelSpacing; info.SliceThickness]';
%                     catch
%                         Spc = [0.4492 0.4492 4]; % reslice it after saving
%                     end
                    count = count+1;
                    temp = dicomread(fullfile(data_dir,sub_dir.name,dir_list(i).name));
                    data.DATA(:,:,count) = squeeze(temp);
                    disp(['Reading image number ' num2str(count) ' of ' num2str(30)])
                end
            end
        end
        
    case 0
        dir_list = dir(data_dir);
        for i = 1:length(dir_list)
            ind(i,:) = ~startsWith(dir_list(i).name,'.') && endsWith(dir_list(i).name,'.dcm');
        end
        
        count = 0;
        for i = 1:length(ind)
            if ind(i) == 1
                if count == 0
                    info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
                    try
                        Spc = [info.PixelSpacing; info.SliceThickness]';
                    catch
                        Spc = [0.5078 0.5078 1]; % reslice it after saving
                    end
                    count = count + 1;
                end
            end
        end
        count = 0;
        for i = 1:length(ind)
            if ind(i) == 1
%                 info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
%                     try
%                         Spc = [info.PixelSpacing; info.SliceThickness]';
%                     catch
%                         Spc = [0.5078 0.5078 1]; % reslice it after saving
%                     end
                count = count+1;
                temp = dicomread(fullfile(data_dir,dir_list(i).name));
                data.DATA(:,:,count) = squeeze(temp);
                disp(['Reading image number ' num2str(count) ' of ' num2str(30)])
            end
        end
        
end

data = data.DATA;
switch save_flag
    case 1
        save(fullfile(tmp.lots_of_dirs{:,1},[date '_T2w_FLAIR.mat']),'data','-mat')
end

end

function [data,Spc] = t1POST(N,directories,save_flag)
disp('T1 POST')
% info = dicominfo('AX_FLAIR/IM-0001-0001-0001.dcm');
% info.ImagePositionPatient;

tmp = directories;
data.DATA = zeros(512,512,170);

temp = dir();
n = size(temp);
for i = 1:n(1)
    if strcmp(fullfile(temp(i).name),'T1w_TFE_postC') ==1
        file_number = i;
    elseif strcmp(fullfile(temp(i).name),'WIP +CsT1W_3D_TFE (WAND)') ==1
        file_number = i;
    elseif strcmp(fullfile(temp(i).name),'WIP AX T1 3D_TFE (WAND)+C 1_2') ==1
        file_number = i;
    elseif strcmp(fullfile(temp(i).name),'WIP AX T1 3D TFE+C') ==1
        file_number = i;
    end
end

data_dir = fullfile(tmp.lots_of_dirs{:,N},fullfile(temp(file_number).name));
sub_dir = dir([data_dir,'/Series*']);

if ~isempty(sub_dir)
    subby = 1;
else
    subby = 0;
end

switch subby
    case 1
        if length(sub_dir) > 1
            sub_dir = sub_dir(end); % this just selects the last directory of the Series because the first directory was used as a test.
        end
        dir_list = dir(fullfile(sub_dir.folder,sub_dir.name));
        for i = 1:length(dir_list)
            % This will open the dicom image that does not have the dcm at
            % the end of it, which is how the data is originally when
            % unzipped. This should handle single files of this name.
            ind(i,:) = ~startsWith(dir_list(i).name,'.') && ~endsWith(dir_list(i).name,'.dcm') && startsWith(dir_list(i).name, 'Image');
        end
        
        if sum(ind) == 0
            for i = 1:length(dir_list)
                % This will open the dicom image that does not have the dcm at
                % the end of it, which is how the data is originally when
                % unzipped. This should handle single files of this name.
                ind(i,:) = ~startsWith(dir_list(i).name,'.') && ~endsWith(dir_list(i).name,'.dcm') && startsWith(dir_list(i).name, 'Mag');
            end
        end
        
        if sum(ind) == 0 % if it still didn't get the files. 
            for i = 1:length(dir_list)
                % This will open the dicom image that does not have the dcm at
                % the end of it, which is how the data is originally when
                % unzipped. This should handle single files of this name.
                ind(i,:) = ~startsWith(dir_list(i).name,'.') && endsWith(dir_list(i).name,'.dcm') && startsWith(dir_list(i).name, 'I');
            end
        end
        
        count = 0;
        for i = 1:length(ind)
            if ind(i) == 1
                if count == 0
                    info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
                    try
                        Spc = [info.PixelSpacing; info.SliceThickness]';
                    catch
                        Spc = [0.5078 0.5078 1]; % reslice it after saving
                    end
                    count = count + 1;
                end
            end
        end
        count = 0;
        if sum(ind) > 1
            temp = dicomread(fullfile(data_dir,sub_dir.name,dir_list(5).name));
            [x,y] = size(temp);
            data.DATA = zeros(x,y,sum(ind));clear temp x y
        end
        if sum(ind) == 1
            for i = 1:length(ind)
                if ind(i) == 1
%                     info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
%                     try
%                         Spc = [info.PixelSpacing; info.SliceThickness]';
%                     catch
%                         Spc = [0.5078 0.5078 1]; % reslice it after saving
%                     end
                    temp = dicomread(fullfile(data_dir,sub_dir.name,dir_list(i).name));
                    data.DATA = squeeze(temp);
                    disp(['Reading image number ' num2str(count) ' of ' num2str(170)])
                end
            end
        else
            for i = 1:length(ind)
                if ind(i) == 1
%                     info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
%                     try
%                         Spc = [info.PixelSpacing; info.SliceThickness]';
%                     catch
%                         Spc = [0.4883 0.4883 1]; % reslice it after saving
%                     end
                    count = count+1;
                    temp = dicomread(fullfile(data_dir,sub_dir.name,dir_list(i).name));
                    data.DATA(:,:,count) = squeeze(temp);
                    disp(['Reading image number ' num2str(count) ' of ' num2str(170)])
                end
            end
        end
        
    case 0
        dir_list = dir(data_dir);
        for i = 1:length(dir_list)
            ind(i,:) = ~startsWith(dir_list(i).name,'.') && endsWith(dir_list(i).name,'.dcm');
        end
        
        count = 0;
        for i = 1:length(ind)
            if ind(i) == 1
                if count == 0
                    info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
                    try
                        Spc = [info.PixelSpacing; info.SliceThickness]';
                    catch
                        Spc = [0.5078 0.5078 1]; % reslice it after saving
                    end
                    count = count + 1;
                end
            end
        end
        count = 0;
        for i = 1:length(ind)
            if ind(i) == 1
%                 info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
%                     try
%                         Spc = [info.PixelSpacing; info.SliceThickness]';
%                     catch
%                         Spc = [0.5078 0.5078 1]; % reslice it after saving
%                     end
                count = count+1;
                temp = dicomread(fullfile(data_dir,dir_list(i).name));
                data.DATA(:,:,count) = squeeze(temp);
                disp(['Reading image number ' num2str(count) ' of ' num2str(170)])
            end
        end
        
end

data = data.DATA;
switch save_flag
    case 1
        save(fullfile(tmp.lots_of_dirs{:,1},[date '_T1w_POST.mat']),'data','-mat')
end

end

function [data,Spc] = SAGE(N,directories,save_flag)
disp('DSC SAGE')
% info = dicominfo('AX_FLAIR/IM-0001-0001-0001.dcm');
% info.ImagePositionPatient;

tmp = directories;
data.DATA = zeros(96,96,11250);

temp = dir();
n = size(temp);
for i = 1:n(1)
    if strcmp(fullfile(temp(i).name),'SAGE_perfusion') ==1
        file_number = i;
    elseif strcmp(fullfile(temp(i).name),'WIP INJECT_SAGEperfusion') ==1
        file_number = i;
    elseif strcmp(fullfile(temp(i).name),'WIP INJECT_MSPERFUSION') ==1
        file_number = i;
    end
end

data_dir = fullfile(tmp.lots_of_dirs{:,N},fullfile(temp(file_number).name));
sub_dir = dir([data_dir,'/Series*']);

if ~isempty(sub_dir)
    subby = 1;
else
    subby = 0;
end

switch subby
    case 1
        if length(sub_dir) > 1
            sub_dir = sub_dir(end); % this just selects the last directory of the Series because the first directory was used as a test.
        end
        dir_list = dir(fullfile(sub_dir.folder,sub_dir.name));
        for i = 1:length(dir_list)
            ind(i,:) = ~startsWith(dir_list(i).name,'.') && ~endsWith(dir_list(i).name,'.dcm') && startsWith(dir_list(i).name, 'Image');
        end
        
        if sum(ind) == 0
            for i = 1:length(dir_list)
                ind(i,:) = ~startsWith(dir_list(i).name,'.') && ~endsWith(dir_list(i).name,'.dcm') && startsWith(dir_list(i).name, 'Mag');
            end
        end
        
        if sum(ind) == 0 % if it still didn't get the files. 
            for i = 1:length(dir_list)
                ind(i,:) = ~startsWith(dir_list(i).name,'.') && endsWith(dir_list(i).name,'.dcm') && startsWith(dir_list(i).name, 'I');
            end
        end
        
        count = 0;
        for i = 1:length(ind)
            if ind(i) == 1
                if count == 0
                    info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
                    try
                        Spc = [info.PixelSpacing; info.SliceThickness]';
                    catch
                        Spc = [0.5078 0.5078 1]; % reslice it after saving
                    end
                    count = count + 1;
                end
            end
        end
        count = 0;
        if sum(ind) > 1
            temp = dicomread(fullfile(data_dir,sub_dir.name,dir_list(5).name));
            [x,y] = size(temp);
            data.DATA = zeros(x,y,sum(ind));clear temp x y
        end
        
        if sum(ind) == 1
            for i = 1:length(ind)
                if ind(i) == 1
                    temp = dicomread(fullfile(data_dir,sub_dir.name,dir_list(i).name));
                    data.DATA = squeeze(temp);
                end
            end
        else
            for i = 1:length(ind)
                if ind(i) == 1
                    count = count+1;
                    temp = dicomread(fullfile(data_dir,sub_dir.name,dir_list(i).name));
                    data.DATA(:,:,count) = squeeze(temp);
                    disp(['Reading image number ' num2str(count) ' of ' num2str(11250)])
                end
            end
        end
        
    case 0
        dir_list = dir(data_dir);
        for i = 1:length(dir_list)
            ind(i,:) = ~startsWith(dir_list(i).name,'.') && endsWith(dir_list(i).name,'.dcm');
        end
        
        count = 0;
        for i = 1:length(ind)
            if ind(i) == 1
                if count == 0
                    info = dicominfo(fullfile(data_dir,sub_dir.name,dir_list(i).name));
                    try
                        Spc = [info.PixelSpacing; info.SliceThickness]';
                    catch
                        Spc = [0.5078 0.5078 1]; % reslice it after saving
                    end
                    count = count + 1;
                end
            end
        end
        count = 0;
        for i = 1:length(ind)
            if ind(i) == 1
                count = count+1;
%                 temp = dicomread(fullfile(data_dir,dir_list(i).name)); % This will not work, files not in order.
                temp = dicomread(fullfile([data_dir,'/IM-0001-' num2str(count,'%04i') '-0001.dcm']));
                data.DATA(:,:,count) = squeeze(temp);
                disp(['Reading image number ' num2str(count) ' of ' num2str(11250)])
            end
        end
        
end

data = data.DATA;
switch save_flag
    case 1
        save(fullfile(tmp.lots_of_dirs{:,1},[date '_SAGE.mat']),'data','-mat')
end

end

function [Data,Parms,AIF,AIFMask,dR2s,Perfusion] = conv_CBV_CBF1(dsc_data_orig)
data = dsc_data_orig;
%    temp = squeeze(reshape(data,[96 96 5 150 15])); Data = permute(temp,[1 2 5 3 4]); % [ 15 5 150 ]
%         temp = squeeze(reshape(data,[96 96 5 15 150])); Data = permute(temp,[1 2 4 3 5]); % [ 15 5 150 ]
%       temp = squeeze(reshape(data,[96 96 15 5 150]));Data = permute(temp,[1 2 3 4 5]); % [ 15 5 150 ]
%       temp = squeeze(reshape(data,[96 96 15 150 5])); Data = permute(temp,[1 2 3 5 4]); % [ 15 5 150 ]
%       temp = squeeze(reshape(data,[96 96 150 15 5])); Data = permute(temp,[1 2 4 5 3]); % [ 15 5 150 ]
      temp = squeeze(reshape(data,[96 96 150 5 15])); Data = permute(temp,[1 2 5 4 3]); % [ 15 5 150 ]
% Data = permute(temp,[1 2 5 3 4]); % [ 15 5 150 ]

Parms.TR = 1800./1000; %TR in s
Parms.time = Parms.TR:Parms.TR:Parms.TR*size(Data,5);
Parms.TE1 = 8.8./1000; %TE in s
Parms.TE2 = 26./1000; %TE in s
Parms.TE5 = 88./1000; %TE in s
Parms.flip= 90; %FA in degrees
[Parms.nx,Parms.ny,Parms.nz,Parms.ne,Parms.nt] = size(Data);
% process data - brain mask and dR2s
Parms.Tes = [8.8 26 nan nan 88]./1000; %NS 20200314
%         filtered_image = abs(cast(DSC.Data,'double'));
filtered_image = (cast(Data,'double'));

[Parms.ss_tp, Parms.prebolus_end_fid,Parms.gd_tp,~] = ...
    determineTimePoints(filtered_image(:,:,:,1:2,:),0,0,Parms.TR); % I'm pretty sure that the vaules for the output are in the wrong order.
prebolus_signal = squeeze(nanmean(filtered_image(:,:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),5)); %average prebolus images together

AIF = AutoAIF_Brain(squeeze(filtered_image(:,:,:,1:2,:)),[Parms.TR ...
    Parms.Tes(1) Parms.Tes(2)],Parms.flip,0,1);title('AIF voxel locations');
load('AIFmask.mat')
AIFMask = AIFmask; % This is not even used.
dR2s.STE0 = squeeze(filtered_image(:,:,:,1,:).*(...
    filtered_image(:,:,:,1,:)./filtered_image(:,:,:,2,:)).^(Parms.TE1/(Parms.TE2-Parms.TE1))); %from Vonken paper
dR2s.STE0_fraction = dR2s.STE0./repmat(nanmean(...
    dR2s.STE0(:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),4),[1 1 1 size(dR2s.STE0,4)]);
% Calculate delta_R2star
prebolus_signal5d = repmat(prebolus_signal,[1 1 1 1 Parms.nt]);
TE_fraction = filtered_image./prebolus_signal5d;
dR2s.all = (1/(Parms.TE2-Parms.TE1)).*squeeze(log(filtered_image(:,:,:,1,:)./...
    filtered_image(:,:,:,2,:)));dR2s.all = dR2s.all - ...
    repmat(nanmean(dR2s.all(:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),4),...
    [1 1 1 size(dR2s.all,4)]);
dR2s.allSE = (1/Parms.TE5)*log(dR2s.STE0./squeeze(filtered_image(:,:,:,5,:)));
dR2s.allSE = dR2s.allSE - repmat(nanmean...
    (dR2s.allSE(:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),4),...
    [1 1 1 size(dR2s.allSE,4)]); %SE DR2 corrected for T1 using STE0!!
%% rCBV
tidx = max(1,round(Parms.gd_tp-60/Parms.TR)):round(Parms.gd_tp+60/Parms.TR);
rtissue = 87;
temp = dR2s.all(:,:,:,tidx)./rtissue;
temp(temp<0) = 0;
Perfusion.rCBV0_all = squeeze(trapz(temp,4))./trapz(AIF(2,tidx));
Perfusion.rCBV0_all(~isfinite(Perfusion.rCBV0_all)) = 0;
temp = dR2s.allSE(:,:,:,tidx)./rtissue;
temp(temp<0) = 0;
Perfusion.rCBV0_allSE = squeeze(trapz(temp,4))./trapz(AIF(2,tidx));
Perfusion.rCBV0_allSE(~isfinite(Perfusion.rCBV0_allSE)) = 0;
%% deconvolution and CBF
dtemp_all = dR2s.all(:,:,:,tidx)./rtissue;
dtemp_allSE = dR2s.allSE(:,:,:,tidx)./rtissue;
zpad = 2;
DR2sAIF = AIF(2,tidx)';
DR2sAIF((length(DR2sAIF)+1):zpad*length(DR2sAIF)) = 0;
if size(DR2sAIF,1)<size(DR2sAIF,2)
    disp('Incorrect - AIF row matrix!');
    DR2sAIF = DR2sAIF';
end
dtemp_all(:,:,:,(size(dtemp_all,4)+1):zpad*size(dtemp_all,4))  = 0;
dtemp_allSE(:,:,:,(size(dtemp_allSE,4)+1):zpad*size(dtemp_allSE,4))  = 0;
AIFmatrixt = Parms.TR*circulant(DR2sAIF,1); %Create circular matrix
%Performs standard SVD on the DeltaR2star values for the AIF
[U,S,V]=svd(AIFmatrixt);
maxS = max(max(S));
S_orig = S;
%Adaptive Threshold for SVD - based on code from Jack Skinner; Liu et al. 1999
Smin = squeeze(nanmin(filtered_image(:,:,:,:,Parms.ss_tp:Parms.prebolus_end_fid+40),[],5)); %find min signal (including during bolus passage)
prebolus_std = squeeze(nanstd(filtered_image(:,:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),[],5)); %find prebolus standard deviation
% % SNRc(j,k)=(S0_TE2/std0_TE2)*(Smin/S0_TE2)*log(S0_TE2/Smin);
SNRc2 = squeeze((prebolus_signal./prebolus_std).*(Smin./prebolus_signal).*log(prebolus_signal./Smin));
SNRc2(isnan(SNRc2)) = 0;
C1 = SNRc2;C2 = SNRc2;C3 = SNRc2;
C1(C1>=4) = 0;C1(C1<2)=0;threshold_C1=(100./(1.7082*SNRc2-1.0988))./100;
C2(C2<4)=0;threshold_C2 = (100./(-0.0061.*(SNRc2.^2)+0.7454*SNRc2+2.8697))./100; %AMS remove >37 criteria! C2(C2>37)=0;
C3(C1~=0)=0;C3(C2~=0)=0;threshold_C3 = 0.15;
threshold = (threshold_C1.*logical(C1))+(threshold_C2.*logical(C2))+(threshold_C3.*logical(C3));

% % %Calculate CBF*R(t)
CBFR_all = zeros(size(dtemp_all));CBFR_allSE = zeros(size(dtemp_allSE));
for i = 1:Parms.nx
    for j = 1:Parms.ny
        for k = 1:Parms.nz
            S = S_orig;
            S(S<(threshold(i,j,k)*maxS))=0;
            %Inverts the diagonal of the S matrix produced through SVD
            S = 1./S;
            S(isinf(S)) = 0;
            %Computes the inverted DeltaR2star values for the AIF
            invDeltaR2starAIF=V*S*U';
            CBFR_all(i,j,k,:) = invDeltaR2starAIF*squeeze(dtemp_all(i,j,k,:));
            CBFR_allSE(i,j,k,:) = invDeltaR2starAIF*squeeze(dtemp_allSE(i,j,k,:));
        end
    end
end

Perfusion.CBF0_all = max(CBFR_all,[],4);
Perfusion.CBF0_all(~isfinite(Perfusion.CBF0_all)) = 0;
Perfusion.rMTT0_all = Perfusion.rCBV0_all./Perfusion.CBF0_all;
Perfusion.rMTT0_all(isinf(Perfusion.rMTT0_all)) = 0;

Perfusion.CBF0_allSE = max(CBFR_allSE,[],4);
Perfusion.CBF0_allSE(~isfinite(Perfusion.CBF0_allSE)) = 0;
Perfusion.rMTT0_allSE = Perfusion.rCBV0_allSE./Perfusion.CBF0_allSE;
Perfusion.rMTT0_allSE(isinf(Perfusion.rMTT0_allSE)) = 0;
end

function [Data,Parms,AIF,AIFMask,dR2s,Perfusion] = conv_CBV_CBF2(dsc_data_orig)
data = dsc_data_orig;
   temp = squeeze(reshape(data,[96 96 5 150 15])); Data = permute(temp,[1 2 5 3 4]); % [ 15 5 150 ]
%         temp = squeeze(reshape(data,[96 96 5 15 150])); Data = permute(temp,[1 2 4 3 5]); % [ 15 5 150 ]
%       temp = squeeze(reshape(data,[96 96 15 5 150]));Data = permute(temp,[1 2 3 4 5]); % [ 15 5 150 ]
%       temp = squeeze(reshape(data,[96 96 15 150 5])); Data = permute(temp,[1 2 3 5 4]); % [ 15 5 150 ]
%       temp = squeeze(reshape(data,[96 96 150 15 5])); Data = permute(temp,[1 2 4 5 3]); % [ 15 5 150 ]
%       temp = squeeze(reshape(data,[96 96 150 5 15])); Data = permute(temp,[1 2 5 4 3]); % [ 15 5 150 ]
% Data = permute(temp,[1 2 5 3 4]); % [ 15 5 150 ]

Parms.TR = 1800./1000; %TR in s
Parms.time = Parms.TR:Parms.TR:Parms.TR*size(Data,5);
Parms.TE1 = 8.8./1000; %TE in s
Parms.TE2 = 26./1000; %TE in s
Parms.TE5 = 88./1000; %TE in s
Parms.flip= 90; %FA in degrees
[Parms.nx,Parms.ny,Parms.nz,Parms.ne,Parms.nt] = size(Data);
% process data - brain mask and dR2s
Parms.Tes = [8.8 26 nan nan 88]./1000; %NS 20200314
%         filtered_image = abs(cast(DSC.Data,'double'));
filtered_image = (cast(Data,'double'));

[Parms.ss_tp, Parms.prebolus_end_fid,Parms.gd_tp,~] = ...
    determineTimePoints(filtered_image(:,:,:,1:2,:),0,0,Parms.TR); % I'm pretty sure that the vaules for the output are in the wrong order.
prebolus_signal = squeeze(nanmean(filtered_image(:,:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),5)); %average prebolus images together

AIF = AutoAIF_Brain(squeeze(filtered_image(:,:,:,1:2,:)),[Parms.TR ...
    Parms.Tes(1) Parms.Tes(2)],Parms.flip,0,1);title('AIF voxel locations');
load('AIFmask.mat')
AIFMask = AIFmask; % This is not even used.
dR2s.STE0 = squeeze(filtered_image(:,:,:,1,:).*(...
    filtered_image(:,:,:,1,:)./filtered_image(:,:,:,2,:)).^(Parms.TE1/(Parms.TE2-Parms.TE1))); %from Vonken paper
dR2s.STE0_fraction = dR2s.STE0./repmat(nanmean(...
    dR2s.STE0(:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),4),[1 1 1 size(dR2s.STE0,4)]);
% Calculate delta_R2star
prebolus_signal5d = repmat(prebolus_signal,[1 1 1 1 Parms.nt]);
TE_fraction = filtered_image./prebolus_signal5d;
dR2s.all = (1/(Parms.TE2-Parms.TE1)).*squeeze(log(filtered_image(:,:,:,1,:)./...
    filtered_image(:,:,:,2,:)));dR2s.all = dR2s.all - ...
    repmat(nanmean(dR2s.all(:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),4),...
    [1 1 1 size(dR2s.all,4)]);
dR2s.allSE = (1/Parms.TE5)*log(dR2s.STE0./squeeze(filtered_image(:,:,:,5,:)));
dR2s.allSE = dR2s.allSE - repmat(nanmean...
    (dR2s.allSE(:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),4),...
    [1 1 1 size(dR2s.allSE,4)]); %SE DR2 corrected for T1 using STE0!!
%% rCBV
tidx = max(1,round(Parms.gd_tp-60/Parms.TR)):round(Parms.gd_tp+60/Parms.TR);
rtissue = 87;
temp = dR2s.all(:,:,:,tidx)./rtissue;
temp(temp<0) = 0;
Perfusion.rCBV0_all = squeeze(trapz(temp,4))./trapz(AIF(2,tidx));
Perfusion.rCBV0_all(~isfinite(Perfusion.rCBV0_all)) = 0;
temp = dR2s.allSE(:,:,:,tidx)./rtissue;
temp(temp<0) = 0;
Perfusion.rCBV0_allSE = squeeze(trapz(temp,4))./trapz(AIF(2,tidx));
Perfusion.rCBV0_allSE(~isfinite(Perfusion.rCBV0_allSE)) = 0;
%% deconvolution and CBF
dtemp_all = dR2s.all(:,:,:,tidx)./rtissue;
dtemp_allSE = dR2s.allSE(:,:,:,tidx)./rtissue;
zpad = 2;
DR2sAIF = AIF(2,tidx)';
DR2sAIF((length(DR2sAIF)+1):zpad*length(DR2sAIF)) = 0;
if size(DR2sAIF,1)<size(DR2sAIF,2)
    disp('Incorrect - AIF row matrix!');
    DR2sAIF = DR2sAIF';
end
dtemp_all(:,:,:,(size(dtemp_all,4)+1):zpad*size(dtemp_all,4))  = 0;
dtemp_allSE(:,:,:,(size(dtemp_allSE,4)+1):zpad*size(dtemp_allSE,4))  = 0;
AIFmatrixt = Parms.TR*circulant(DR2sAIF,1); %Create circular matrix
%Performs standard SVD on the DeltaR2star values for the AIF
[U,S,V]=svd(AIFmatrixt);
maxS = max(max(S));
S_orig = S;
%Adaptive Threshold for SVD - based on code from Jack Skinner; Liu et al. 1999
Smin = squeeze(nanmin(filtered_image(:,:,:,:,Parms.ss_tp:Parms.prebolus_end_fid+40),[],5)); %find min signal (including during bolus passage)
prebolus_std = squeeze(nanstd(filtered_image(:,:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),[],5)); %find prebolus standard deviation
% % SNRc(j,k)=(S0_TE2/std0_TE2)*(Smin/S0_TE2)*log(S0_TE2/Smin);
SNRc2 = squeeze((prebolus_signal./prebolus_std).*(Smin./prebolus_signal).*log(prebolus_signal./Smin));
SNRc2(isnan(SNRc2)) = 0;
C1 = SNRc2;C2 = SNRc2;C3 = SNRc2;
C1(C1>=4) = 0;C1(C1<2)=0;threshold_C1=(100./(1.7082*SNRc2-1.0988))./100;
C2(C2<4)=0;threshold_C2 = (100./(-0.0061.*(SNRc2.^2)+0.7454*SNRc2+2.8697))./100; %AMS remove >37 criteria! C2(C2>37)=0;
C3(C1~=0)=0;C3(C2~=0)=0;threshold_C3 = 0.15;
threshold = (threshold_C1.*logical(C1))+(threshold_C2.*logical(C2))+(threshold_C3.*logical(C3));

% % %Calculate CBF*R(t)
CBFR_all = zeros(size(dtemp_all));CBFR_allSE = zeros(size(dtemp_allSE));
for i = 1:Parms.nx
    for j = 1:Parms.ny
        for k = 1:Parms.nz
            S = S_orig;
            S(S<(threshold(i,j,k)*maxS))=0;
            %Inverts the diagonal of the S matrix produced through SVD
            S = 1./S;
            S(isinf(S)) = 0;
            %Computes the inverted DeltaR2star values for the AIF
            invDeltaR2starAIF=V*S*U';
            CBFR_all(i,j,k,:) = invDeltaR2starAIF*squeeze(dtemp_all(i,j,k,:));
            CBFR_allSE(i,j,k,:) = invDeltaR2starAIF*squeeze(dtemp_allSE(i,j,k,:));
        end
    end
end

Perfusion.CBF0_all = max(CBFR_all,[],4);
Perfusion.CBF0_all(~isfinite(Perfusion.CBF0_all)) = 0;
Perfusion.rMTT0_all = Perfusion.rCBV0_all./Perfusion.CBF0_all;
Perfusion.rMTT0_all(isinf(Perfusion.rMTT0_all)) = 0;

Perfusion.CBF0_allSE = max(CBFR_allSE,[],4);
Perfusion.CBF0_allSE(~isfinite(Perfusion.CBF0_allSE)) = 0;
Perfusion.rMTT0_allSE = Perfusion.rCBV0_allSE./Perfusion.CBF0_allSE;
Perfusion.rMTT0_allSE(isinf(Perfusion.rMTT0_allSE)) = 0;
end

function [Data,Parms,AIF,AIFMask,dR2s,Perfusion] = conv_CBV_CBF3(dsc_data_orig)
data = dsc_data_orig;
        %         data = imrotate3(data,-90,[0 0 1]);
        %         temp = squeeze(reshape(data,[96 96 15 150 5])); Data = permute(temp,[1 2 3 5 4]); % [ 15 5 150 ]
        temp = squeeze(reshape(data,[96 96 5 150 15])); Data = permute(temp,[1 2 5 3 4]); % [ 15 5 150 ] %this one worked for loading from single Image (0001)
        %         temp = squeeze(reshape(data,[96 96 15 5 150]));Data = permute(temp,[1 2 3 4 5]); % [ 15 5 150 ]
        %         temp = squeeze(reshape(data,[96 96 150 15 5])); Data = permute(temp,[1 2 5 3 4]); % [ 15 5 150 ]
        %         temp = squeeze(reshape(data,[96 96 150 5 15])); Data = permute(temp,[1 2 5 4 3]); % [ 15 5 150 ] % This one worked if the data is from a file of many dicoms

Parms.TR = 1800./1000; %TR in s
Parms.time = Parms.TR:Parms.TR:Parms.TR*size(Data,5);
Parms.TE1 = 8.8./1000; %TE in s
Parms.TE2 = 26./1000; %TE in s
Parms.TE5 = 88./1000; %TE in s
Parms.flip= 90; %FA in degrees
[Parms.nx,Parms.ny,Parms.nz,Parms.ne,Parms.nt] = size(Data);
% process data - brain mask and dR2s
Parms.Tes = [8.8 26 nan nan 88]./1000; %NS 20200314
%         filtered_image = abs(cast(DSC.Data,'double'));
filtered_image = (cast(Data,'double'));

[Parms.ss_tp, Parms.prebolus_end_fid,Parms.gd_tp,~] = ...
    determineTimePoints(filtered_image(:,:,:,1:2,:),0,0,Parms.TR); % I'm pretty sure that the vaules for the output are in the wrong order.
prebolus_signal = squeeze(nanmean(filtered_image(:,:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),5)); %average prebolus images together

AIF = AutoAIF_Brain(squeeze(filtered_image(:,:,:,1:2,:)),[Parms.TR ...
    Parms.Tes(1) Parms.Tes(2)],Parms.flip,0,1);title('AIF voxel locations');
load('AIFmask.mat')
AIFMask = AIFmask; % This is not even used.
dR2s.STE0 = squeeze(filtered_image(:,:,:,1,:).*(...
    filtered_image(:,:,:,1,:)./filtered_image(:,:,:,2,:)).^(Parms.TE1/(Parms.TE2-Parms.TE1))); %from Vonken paper
dR2s.STE0_fraction = dR2s.STE0./repmat(nanmean(...
    dR2s.STE0(:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),4),[1 1 1 size(dR2s.STE0,4)]);
% Calculate delta_R2star
prebolus_signal5d = repmat(prebolus_signal,[1 1 1 1 Parms.nt]);
TE_fraction = filtered_image./prebolus_signal5d;
dR2s.all = (1/(Parms.TE2-Parms.TE1)).*squeeze(log(filtered_image(:,:,:,1,:)./...
    filtered_image(:,:,:,2,:)));dR2s.all = dR2s.all - ...
    repmat(nanmean(dR2s.all(:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),4),...
    [1 1 1 size(dR2s.all,4)]);
dR2s.allSE = (1/Parms.TE5)*log(dR2s.STE0./squeeze(filtered_image(:,:,:,5,:)));
dR2s.allSE = dR2s.allSE - repmat(nanmean...
    (dR2s.allSE(:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),4),...
    [1 1 1 size(dR2s.allSE,4)]); %SE DR2 corrected for T1 using STE0!!
%% rCBV
tidx = max(1,round(Parms.gd_tp-60/Parms.TR)):round(Parms.gd_tp+60/Parms.TR);
rtissue = 87;
temp = dR2s.all(:,:,:,tidx)./rtissue;
temp(temp<0) = 0;
Perfusion.rCBV0_all = squeeze(trapz(temp,4))./trapz(AIF(2,tidx));
Perfusion.rCBV0_all(~isfinite(Perfusion.rCBV0_all)) = 0;
temp = dR2s.allSE(:,:,:,tidx)./rtissue;
temp(temp<0) = 0;
Perfusion.rCBV0_allSE = squeeze(trapz(temp,4))./trapz(AIF(2,tidx));
Perfusion.rCBV0_allSE(~isfinite(Perfusion.rCBV0_allSE)) = 0;
%% deconvolution and CBF
dtemp_all = dR2s.all(:,:,:,tidx)./rtissue;
dtemp_allSE = dR2s.allSE(:,:,:,tidx)./rtissue;
zpad = 2;
DR2sAIF = AIF(2,tidx)';
DR2sAIF((length(DR2sAIF)+1):zpad*length(DR2sAIF)) = 0;
if size(DR2sAIF,1)<size(DR2sAIF,2)
    disp('Incorrect - AIF row matrix!');
    DR2sAIF = DR2sAIF';
end
dtemp_all(:,:,:,(size(dtemp_all,4)+1):zpad*size(dtemp_all,4))  = 0;
dtemp_allSE(:,:,:,(size(dtemp_allSE,4)+1):zpad*size(dtemp_allSE,4))  = 0;
AIFmatrixt = Parms.TR*circulant(DR2sAIF,1); %Create circular matrix
%Performs standard SVD on the DeltaR2star values for the AIF
[U,S,V]=svd(AIFmatrixt);
maxS = max(max(S));
S_orig = S;
%Adaptive Threshold for SVD - based on code from Jack Skinner; Liu et al. 1999
Smin = squeeze(nanmin(filtered_image(:,:,:,:,Parms.ss_tp:Parms.prebolus_end_fid+40),[],5)); %find min signal (including during bolus passage)
prebolus_std = squeeze(nanstd(filtered_image(:,:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),[],5)); %find prebolus standard deviation
% % SNRc(j,k)=(S0_TE2/std0_TE2)*(Smin/S0_TE2)*log(S0_TE2/Smin);
SNRc2 = squeeze((prebolus_signal./prebolus_std).*(Smin./prebolus_signal).*log(prebolus_signal./Smin));
SNRc2(isnan(SNRc2)) = 0;
C1 = SNRc2;C2 = SNRc2;C3 = SNRc2;
C1(C1>=4) = 0;C1(C1<2)=0;threshold_C1=(100./(1.7082*SNRc2-1.0988))./100;
C2(C2<4)=0;threshold_C2 = (100./(-0.0061.*(SNRc2.^2)+0.7454*SNRc2+2.8697))./100; %AMS remove >37 criteria! C2(C2>37)=0;
C3(C1~=0)=0;C3(C2~=0)=0;threshold_C3 = 0.15;
threshold = (threshold_C1.*logical(C1))+(threshold_C2.*logical(C2))+(threshold_C3.*logical(C3));

% % %Calculate CBF*R(t)
CBFR_all = zeros(size(dtemp_all));CBFR_allSE = zeros(size(dtemp_allSE));
for i = 1:Parms.nx
    for j = 1:Parms.ny
        for k = 1:Parms.nz
            S = S_orig;
            S(S<(threshold(i,j,k)*maxS))=0;
            %Inverts the diagonal of the S matrix produced through SVD
            S = 1./S;
            S(isinf(S)) = 0;
            %Computes the inverted DeltaR2star values for the AIF
            invDeltaR2starAIF=V*S*U';
            CBFR_all(i,j,k,:) = invDeltaR2starAIF*squeeze(dtemp_all(i,j,k,:));
            CBFR_allSE(i,j,k,:) = invDeltaR2starAIF*squeeze(dtemp_allSE(i,j,k,:));
        end
    end
end

Perfusion.CBF0_all = max(CBFR_all,[],4);
Perfusion.CBF0_all(~isfinite(Perfusion.CBF0_all)) = 0;
Perfusion.rMTT0_all = Perfusion.rCBV0_all./Perfusion.CBF0_all;
Perfusion.rMTT0_all(isinf(Perfusion.rMTT0_all)) = 0;

Perfusion.CBF0_allSE = max(CBFR_allSE,[],4);
Perfusion.CBF0_allSE(~isfinite(Perfusion.CBF0_allSE)) = 0;
Perfusion.rMTT0_allSE = Perfusion.rCBV0_allSE./Perfusion.CBF0_allSE;
Perfusion.rMTT0_allSE(isinf(Perfusion.rMTT0_allSE)) = 0;
end

function [Data,Parms,AIF,AIFMask,dR2s,Perfusion] = conv_CBV_CBF4(dsc_data_orig)
data = dsc_data_orig;
        %         data = imrotate3(data,-90,[0 0 1]);
        %         temp = squeeze(reshape(data,[96 96 15 150 5])); Data = permute(temp,[1 2 3 5 4]); % [ 15 5 150 ]
%         temp = squeeze(reshape(data,[96 96 5 150 15])); Data = permute(temp,[1 2 5 3 4]); % [ 15 5 150 ] %this one worked for loading from single Image (0001)
        %         temp = squeeze(reshape(data,[96 96 15 5 150]));Data = permute(temp,[1 2 3 4 5]); % [ 15 5 150 ]
        %         temp = squeeze(reshape(data,[96 96 150 15 5])); Data = permute(temp,[1 2 5 3 4]); % [ 15 5 150 ]
                temp = squeeze(reshape(data,[96 96 150 5 15])); Data = permute(temp,[1 2 5 4 3]); % [ 15 5 150 ] % This one worked if the data is from a file of many dicoms

Parms.TR = 1800./1000; %TR in s
Parms.time = Parms.TR:Parms.TR:Parms.TR*size(Data,5);
Parms.TE1 = 8.8./1000; %TE in s
Parms.TE2 = 26./1000; %TE in s
Parms.TE5 = 88./1000; %TE in s
Parms.flip= 90; %FA in degrees
[Parms.nx,Parms.ny,Parms.nz,Parms.ne,Parms.nt] = size(Data);
% process data - brain mask and dR2s
Parms.Tes = [8.8 26 nan nan 88]./1000; %NS 20200314
%         filtered_image = abs(cast(DSC.Data,'double'));
filtered_image = (cast(Data,'double'));

[Parms.ss_tp, Parms.prebolus_end_fid,Parms.gd_tp,~] = ...
    determineTimePoints(filtered_image(:,:,:,1:2,:),0,0,Parms.TR); % I'm pretty sure that the vaules for the output are in the wrong order.
prebolus_signal = squeeze(nanmean(filtered_image(:,:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),5)); %average prebolus images together

AIF = AutoAIF_Brain(squeeze(filtered_image(:,:,:,1:2,:)),[Parms.TR ...
    Parms.Tes(1) Parms.Tes(2)],Parms.flip,0,1);title('AIF voxel locations');
load('AIFmask.mat')
AIFMask = AIFmask; % This is not even used.
dR2s.STE0 = squeeze(filtered_image(:,:,:,1,:).*(...
    filtered_image(:,:,:,1,:)./filtered_image(:,:,:,2,:)).^(Parms.TE1/(Parms.TE2-Parms.TE1))); %from Vonken paper
dR2s.STE0_fraction = dR2s.STE0./repmat(nanmean(...
    dR2s.STE0(:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),4),[1 1 1 size(dR2s.STE0,4)]);
% Calculate delta_R2star
prebolus_signal5d = repmat(prebolus_signal,[1 1 1 1 Parms.nt]);
TE_fraction = filtered_image./prebolus_signal5d;
dR2s.all = (1/(Parms.TE2-Parms.TE1)).*squeeze(log(filtered_image(:,:,:,1,:)./...
    filtered_image(:,:,:,2,:)));dR2s.all = dR2s.all - ...
    repmat(nanmean(dR2s.all(:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),4),...
    [1 1 1 size(dR2s.all,4)]);
dR2s.allSE = (1/Parms.TE5)*log(dR2s.STE0./squeeze(filtered_image(:,:,:,5,:)));
dR2s.allSE = dR2s.allSE - repmat(nanmean...
    (dR2s.allSE(:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),4),...
    [1 1 1 size(dR2s.allSE,4)]); %SE DR2 corrected for T1 using STE0!!
%% rCBV
tidx = max(1,round(Parms.gd_tp-60/Parms.TR)):round(Parms.gd_tp+60/Parms.TR);
rtissue = 87;
temp = dR2s.all(:,:,:,tidx)./rtissue;
temp(temp<0) = 0;
Perfusion.rCBV0_all = squeeze(trapz(temp,4))./trapz(AIF(2,tidx));
Perfusion.rCBV0_all(~isfinite(Perfusion.rCBV0_all)) = 0;
temp = dR2s.allSE(:,:,:,tidx)./rtissue;
temp(temp<0) = 0;
Perfusion.rCBV0_allSE = squeeze(trapz(temp,4))./trapz(AIF(2,tidx));
Perfusion.rCBV0_allSE(~isfinite(Perfusion.rCBV0_allSE)) = 0;
%% deconvolution and CBF
dtemp_all = dR2s.all(:,:,:,tidx)./rtissue;
dtemp_allSE = dR2s.allSE(:,:,:,tidx)./rtissue;
zpad = 2;
DR2sAIF = AIF(2,tidx)';
DR2sAIF((length(DR2sAIF)+1):zpad*length(DR2sAIF)) = 0;
if size(DR2sAIF,1)<size(DR2sAIF,2)
    disp('Incorrect - AIF row matrix!');
    DR2sAIF = DR2sAIF';
end
dtemp_all(:,:,:,(size(dtemp_all,4)+1):zpad*size(dtemp_all,4))  = 0;
dtemp_allSE(:,:,:,(size(dtemp_allSE,4)+1):zpad*size(dtemp_allSE,4))  = 0;
AIFmatrixt = Parms.TR*circulant(DR2sAIF,1); %Create circular matrix
%Performs standard SVD on the DeltaR2star values for the AIF
[U,S,V]=svd(AIFmatrixt);
maxS = max(max(S));
S_orig = S;
%Adaptive Threshold for SVD - based on code from Jack Skinner; Liu et al. 1999
Smin = squeeze(nanmin(filtered_image(:,:,:,:,Parms.ss_tp:Parms.prebolus_end_fid+40),[],5)); %find min signal (including during bolus passage)
prebolus_std = squeeze(nanstd(filtered_image(:,:,:,:,Parms.ss_tp:Parms.prebolus_end_fid),[],5)); %find prebolus standard deviation
% % SNRc(j,k)=(S0_TE2/std0_TE2)*(Smin/S0_TE2)*log(S0_TE2/Smin);
SNRc2 = squeeze((prebolus_signal./prebolus_std).*(Smin./prebolus_signal).*log(prebolus_signal./Smin));
SNRc2(isnan(SNRc2)) = 0;
C1 = SNRc2;C2 = SNRc2;C3 = SNRc2;
C1(C1>=4) = 0;C1(C1<2)=0;threshold_C1=(100./(1.7082*SNRc2-1.0988))./100;
C2(C2<4)=0;threshold_C2 = (100./(-0.0061.*(SNRc2.^2)+0.7454*SNRc2+2.8697))./100; %AMS remove >37 criteria! C2(C2>37)=0;
C3(C1~=0)=0;C3(C2~=0)=0;threshold_C3 = 0.15;
threshold = (threshold_C1.*logical(C1))+(threshold_C2.*logical(C2))+(threshold_C3.*logical(C3));

% % %Calculate CBF*R(t)
CBFR_all = zeros(size(dtemp_all));CBFR_allSE = zeros(size(dtemp_allSE));
for i = 1:Parms.nx
    for j = 1:Parms.ny
        for k = 1:Parms.nz
            S = S_orig;
            S(S<(threshold(i,j,k)*maxS))=0;
            %Inverts the diagonal of the S matrix produced through SVD
            S = 1./S;
            S(isinf(S)) = 0;
            %Computes the inverted DeltaR2star values for the AIF
            invDeltaR2starAIF=V*S*U';
            CBFR_all(i,j,k,:) = invDeltaR2starAIF*squeeze(dtemp_all(i,j,k,:));
            CBFR_allSE(i,j,k,:) = invDeltaR2starAIF*squeeze(dtemp_allSE(i,j,k,:));
        end
    end
end

Perfusion.CBF0_all = max(CBFR_all,[],4);
Perfusion.CBF0_all(~isfinite(Perfusion.CBF0_all)) = 0;
Perfusion.rMTT0_all = Perfusion.rCBV0_all./Perfusion.CBF0_all;
Perfusion.rMTT0_all(isinf(Perfusion.rMTT0_all)) = 0;

Perfusion.CBF0_allSE = max(CBFR_allSE,[],4);
Perfusion.CBF0_allSE(~isfinite(Perfusion.CBF0_allSE)) = 0;
Perfusion.rMTT0_allSE = Perfusion.rCBV0_allSE./Perfusion.CBF0_allSE;
Perfusion.rMTT0_allSE(isinf(Perfusion.rMTT0_allSE)) = 0;
end

function [none] = run_flirt(in_file,out_file,ref_file,omat_file)

fsl = getenv('FSLDIR'); %set the path for FSL
dirname = pwd;

in = in_file;
out = out_file;
ref = ref_file;
omat = omat_file;

disp(['Registration using ' in ' ' out ' ' ref ' ' omat])

opts = ' -bins 256 -cost normmi -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear ';
opts_files = [' -in ' in ' -ref ' ref ' -out ' out ' -omat ' omat];


unix([fsl '/bin/flirt' opts_files opts]);

disp('Registration Finished')

end
