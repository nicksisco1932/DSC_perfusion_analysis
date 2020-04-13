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
    save_file_flag = 1; % easier than typing it for every function.
%     try
%         [T2w_FLAIR.Data, T2w_FLAIR.Spc] = t2FLAIR(i,tmp,save_file_flag );
%         tempI = make_nii(imrotate3(T2w_FLAIR.Data,-90,[0 0 1]),T2w_FLAIR.Spc,[0 0 0]);
%         save_nii(tempI, [pwd  '/T2w_FLAIR.nii'], 0);
%     catch
%         warning(['Check on ' a{end} 'T2w_FLAIR']);
%     end
%     try
%         [T1w_pre.Data, T1w_pre.Spc] = t1PRE(i,tmp,save_file_flag );
%         tempI = make_nii(imrotate3(T1w_pre.Data,-90,[0 0 1]),T1w_pre.Spc,[0 0 0]);
%         save_nii(tempI, [pwd  '/T1w_pre.nii'], 0);
%     catch
%         warning(['Check on ' a{end} ' T1w_PREC']);
%     end
%     try
%         [T1w_post.Data, T1w_post.Spc] = t1POST(i,tmp,save_file_flag );
%         tempI = make_nii(imrotate3(T1w_post.Data,-90,[0 0 1]),T1w_post.Spc,[0 0 0]);
%         save_nii(tempI, [pwd  '/T1w_post.nii'], 0);
%     catch
%         warning(['Check on ' a{end} 'T1w_POSTC']);
%     end
    try
        [DSC.Data, DSC.Spc] = SAGE(i,tmp,save_file_flag );
%         tempI = make_nii(imrotate3(DSC.Data,-90,[0 0 1]),DSC.Spc,[0 0 0]);
%         save_nii(tempI, [pwd  '/SAGE.nii'], 0);
    catch
        warning(['Check on ' a{end} 'DSC']);
    end

   
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
                        Spc = [2.5 2.5 5]; % reslice it after saving
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
                    temp = dicomread(fullfile([data_dir,'/' sub_dir.name,'/Image (' num2str(count,'%04i') ')' ]));
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
                        Spc = [2.5 2.5 5]; % reslice it after saving
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
        save(fullfile(tmp.lots_of_dirs{:,N},[date '_SAGE.mat']),'data','-mat')
    case 0
        disp('not saving')
end

end
