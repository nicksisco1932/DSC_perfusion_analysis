%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                             General Perfusion Parameters                                                 %
% This function determines steady state and contrast arrival time points for a perfusion volume, and masks the brain.      %
%                                                                                                                          %
% Laura Bell 1/5/2016                                                                                                      %
%                                                                                                                          %
% Usage: [ss_tp, gd_tp, pk_tp] = determineTimePoints(vol, VOXflag, SaveFlag, TR);                                          %
% Inputs:                                                                                                                  %
% - vol is a 4D [nx ny nz nt] or 5D [nx ny nz ne nt] matrix                                                                %
% - VOXflag = 0; determine time points on average brain signal                                                             %
%           = 1; determine time points on a voxel by voxel basis                                                           %
% - SaveFlag = 0; don't save outputs                                                                                       %
%            = 1; save outputs                                                                                             %
%                                                                                                                          %
% Outputs:                                                                                                                 %
% - ss_tp = time point at which steady state has been reached (if dummy scans are included)                                %
% - gd_tp = time point at which gadolium arrives                                                                           %
% - pk_tp = time point at peak gadolium concentration                                                                      %
%                                                                                                                          %
% Additional functions needed to run this function:                                                                        %
% - Region Growing by Christian Wuerslin                                                                                   %
% http://www.mathworks.com/matlabcentral/fileexchange/41666-fast-3d-2d-region-growing--mex-                                %
%                                                                                                                          %
% Edits:                                                                                                                   %
% - LCB 2/5/2016                                                                                                           %
%       Changed assignment of SaveFlag (=1 should mean save!)                                                              %
%       Changed smooth option to avoid bright signal due to fluid in resection volumes                                     %
%       Added imfill function to fill in holes within the brain                                                            %
%       Added option to output brainmask to workspace                                                                      %
% - LCB 4/6/2016                                                                                                           %
%       Added option to output mean of signal for the second echo if exists                                                %
%       Changed location of timepoint arrows for figure                                                                    %
% - LCB 7/1/2016                                                                                                           %
%       Dilutes (imdilate) the brainmask as a final step to remove erroneous fat around the skull                          %
%       In ss_tp determined changed slope >= 0 to -1 (works better for single echo EPI)                                    %
%       Added TR input as an option, otherwise assumes 1 sec                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ss_tp, gd_tp, pk_tp, brainmask] = determineTimePoints(vol, VOXflag, SaveFlag, TR)
if length(size(vol)) == 5
    echoFlag = 1;
elseif length(size(vol)) == 4
    echoFlag = 0;
else
    warning('Incorrect Matrix Size, this might be a time series.') %NS 20200324
    return;
end


%LCB 7/1/2016
if ~exist('TR','var')
    TR = 1;
    warning('assuming TR is 1 second');
end
%% First mask brain and determine steady state and contrast arrival timepoints
switch echoFlag
    case 1
        [nx,ny,nz,ne,nt] = size(vol);
        I = squeeze(vol(:,:,:,2,:)); %only need first echo dataset for this section
        I_v = reshape(I, [nx*ny*nz nt]);
    otherwise
        [nx,ny,nz,nt] = size(vol);
        ne = 1; %LCB 4/6/2016
        I = vol;
        I_v = reshape(I, [nx*ny*nz nt]);
end

% Mask brain by signal threshold and region growing function
maskmx = max(I_v(:,1));
I_vnorm = I_v(:,1)/maskmx*100;

h(1) = figure('Color', [1 1 1],'Visible', 'on');
hs = histogram(I_vnorm);
title('Histogram of Image SI normalized to max SI');
hold on;
xlabel('Normalized SI', 'FontSize', 14);
ylabel('Counts', 'FontSize', 14)
xlimits = hs.BinLimits; binwidth = hs.BinWidth;
x_axis = xlimits(1)+binwidth:binwidth:xlimits(2);

[pks,locs,w,p] = findpeaks(smooth(hs.Values,'rloess'), x_axis);
table = cat(2, locs', w', p, pks);
[~,mx_w] = max(w); %assume biggest width location
threshold = table(mx_w,1) - table(mx_w,2);
I_vnorm(I_vnorm < threshold) = 0;
I_vnorm(I_vnorm > 0) = 1; %make binary
masktmp = reshape(I_vnorm, [nx ny nz]);
ylim = get(gca, 'ylim');
line([floor(threshold) floor(threshold)], ylim, 'Color', 'r', 'LineWidth', 2)
str1 = '\rightarrow Keep these pixels';
text(floor(threshold), max(ylim)*0.5, str1,'FontSize',14);

brainmask = zeros(size(masktmp));
for i = 1:nz
    slice = masktmp(:,:,i);
    CC = bwconncomp(slice);
    
    %determine number of regions - max mostly 2 when brain hemispheres
    %split on cephalic slices
    numPixels = cellfun(@numel,CC.PixelIdxList);
    locObjects = find(numPixels./max(numPixels) > 0.9);
    numObjects = length(locObjects);
    tmpbrainmask = zeros(nx,ny,numObjects);
    for indexObject = 1:numObjects
        slice(CC.PixelIdxList{locObjects(indexObject)}) = 2;
        slice(slice ~= 2) = 0; slice(slice == 2) = 1; %make binary
        COMstruct = regionprops(slice, squeeze(I(:,:,i,1)), {'Centroid','WeightedCentroid'});
        COM = cat(1, COMstruct.WeightedCentroid);
        filledslice = imfill(slice);
        if filledslice(round(COM(2)), round(COM(1))) == 1
            tmpbrainmask(:,:,indexObject) = RegionGrowing(imfill(slice), 1, [round(COM(2)) round(COM(1))]); %NS change 20200220
        elseif filledslice(ceil(COM(2)), ceil(COM(1))) == 1
            tmpbrainmask(:,:,indexObject) = RegionGrowing(imfill(slice), 1, [ceil(COM(2)) ceil(COM(1))]);
        elseif filledslice(floor(COM(2)), floor(COM(1))) == 1
            tmpbrainmask(:,:,indexObject) = RegionGrowing(imfill(slice), 1, [floor(COM(2)) floor(COM(1))]);
        end
        %         tmpbrainmask(:,:,indexObject)=RegionGrowing(imfill(slice),1,[round(COM(2)) round(COM(1))])
    end
    brainmask(:,:,i) = sum(tmpbrainmask,3);
end


%LCB 7/1/2016
er_brainmask = zeros(nx,ny,nz);
for z = 1:nz
    er_brainmask(:,:,z) = imerode(imfill(brainmask(:,:,z)),strel('disk', 2, 0)); %AMS changed strel('disk',5,0) to strel('disk',2,0), bigger mask
end
%figure, montage(permute(er_brainmask, [1 2 4 3])); title('Brain mask');
brainmask = er_brainmask;

if isequal(SaveFlag, 1)
    save brainmask.mat brainmask
end

h(2) = figure('Color', [1 1 1],'Visible', 'on');
montage(reshape(brainmask, [nx ny 1 nz]));
text(20, 20, 'Final Brain Mask By Slice', 'Color', 'w', 'FontSize', 14);
% clear masktmp xx COM COMloc slicemask maskmx masktmp hs % don't have to clear things in a function

switch VOXflag
    case 0
        %Find steady state and contrast arrival timepoints
        mean_tc = mean(reshape(I, [nx*ny*nz nt])); %mean time course. This reshape is not neccessary since MATLAB indeces would take care of it
        slope = floor(diff(mean_tc)); %find derivative approximate of curve to first determine steady state location
        % %         if slope(1) < -20 %AMS commented out 7/21/2016
        % %             ss_tp = find(slope >= -1, 1, 'first') + round(3/TR); %steady state reached, add 3 seconds for good measure; LCB 7/1/2016 changed slope >= 0 to -1 %AMS commented out 7/21/2016
        if slope(1) < -20 || slope(1) > 20 %include positive slope - important for SAGE
            ss_tp = find(slope >= -1 & slope <= 1, 1, 'first') + round(2/TR); %steady state reached, add 2 seconds for good measure; LCB 7/1/2016 changed slope >= 0 to -1; %AMS 7/21/2016 -1 <= slope <= 1
            if ss_tp > nt/2
                warning('Did not find steady-state in first half of time-course, trying again')
                ss_tp = find(slope >= -20 & slope <= 20, 1, 'first') + round(2/TR); %steady state reached, add 2 seconds for good measure; LCB 7/1/2016 changed slope >= 0 to -1; %AMS 7/21/2016 -1 <= slope <= 1
            end
        else
            ss_tp = 1; %assume no dummy scans included in data
        end
        if isempty(ss_tp) %catch empty ss_tp??
            ss_tp = 1;
        end
        %         [pks,locs,w] = findpeaks(-mean_tc(:,ss_tp:end));
        [pks,locs,w] = findpeaks(-smoothdata(mean_tc(:,ss_tp:end)));
        [~, gd_index] = max(pks);
        woof = diff(mean_tc);
        gd_index = find(abs(woof) == max(abs(woof)));
        
        woof = (diff(diff(-mean_tc)));
        gd_index = find(woof ==max(diff(diff(-mean_tc))));
        gd_tp = gd_index;
%         gd_tp = locs(gd_index) + ss_tp + 1 - round(3/TR) - floor(w(gd_index)*1.5); %contrast arrival % changing this NS
        
%         pk_tp = locs(gd_index) + ss_tp - round(3/TR) + 2; %time to peak
        pk_tp = find(mean_tc == min(mean_tc));
        h(3) = figure('Color', [1 1 1]);
        title('Mean Time Course Across Volume'); hold on;
        plot(mean_tc, 'k', 'LineWidth', 2'); hold on;
        
        if ne == 2 %LCB 4/6/2016 % ne is the number of echos NS
            I2 = squeeze(vol(:,:,:,1,:)); %only need first echo dataset for this section
            mean_tc2 = mean(reshape(I2, [nx*ny*nz nt])); %mean time course
            plot(mean_tc2, '--k', 'LineWidth', 2'); hold on;
            legend('Second Echo', 'First Echo'); %NS changed, pretty sure this was reversed on accident.
            %         elseif ne == 1 % NS 20200402, might not work
            %             I2 = squeeze(vol(:,:,:,1,:)); %only need first echo dataset for this section
            %             mean_tc2 = mean(reshape(I2, [nx*ny*nz nt])); %mean time course
            %             plot(mean_tc2, '--k', 'LineWidth', 2'); hold on;
            %             legend('First Echo');
        end
        
        xlabel('Time points', 'FontSize', 14); ylabel('Signal A.U.', 'FontSize', 14);
        ylim = get(gca, 'ylim');
        line([ss_tp ss_tp], ylim, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
        line([gd_tp gd_tp], ylim, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
        str1 = '\rightarrow Steady State Reached';
        text(ss_tp, mean_tc(ss_tp)*1.25, str1, 'FontSize',14); %LCB 4/6/2016
        str1 = '\rightarrow Contrast Arrival';
        text(gd_tp, mean_tc(gd_tp)*1.1, str1,'FontSize',14); %LCB 4/6/2016
        
        if isequal(SaveFlag, 1)
            savefig(h, 'Timepoint_steps.fig');
        end
        
    case 1
        warning('Needs to be coded.'); % What? NS 20200324 I'm changing gears on this a bit
%         ind = find(brainmask == 1);
%         temp = reshape(I,[nx*ny*nz nt]);
%         integrated = trapz(temp)./max(trapz(temp)); %normalized to the max signal
% %         for i = 1:ind % I think this for loop needs to feed into this
% %         script since the values are not saved inside a function.
%         for i = 483
%             %             plot(1:150,temp(ind(i),:))
%             tempY = temp(ind(i),:); % This will need to change.
%             %             slope = floor(diff(temp(ind(i),:)));
%             slope = diff(smoothdata(tempY,'gaussian'));
% %             if slope(1) < -20 || slope(1) > 20 %include positive slope - important for SAGE
% %                 ss_tp = find(slope >= -1 & slope <= 1, 1, 'first') + round(2/TR); %steady state reached, add 2 seconds for good measure; LCB 7/1/2016 changed slope >= 0 to -1; %AMS 7/21/2016 -1 <= slope <= 1
% %                 if ss_tp > nt/2
% %                     warning('Did not find steady-state in first half of time-course, trying again')
% %                     ss_tp = find(slope >= -20 & slope <= 20, 1, 'first') + round(2/TR); %steady state reached, add 2 seconds for good measure; LCB 7/1/2016 changed slope >= 0 to -1; %AMS 7/21/2016 -1 <= slope <= 1
% %                 end
% %             else
% %                 ss_tp = 1; %assume no dummy scans included in data
% %             end
% %             if isempty(ss_tp) %catch empty ss_tp??
% %                 ss_tp = 1;
% %             end
%             ss_tp = 1; % just going to hard code this here for now since the above code is a mess and it mostly differs to 1
%             [pks,locs,w] = findpeaks(-tempY(:,ss_tp:end));
%             i
% %             [~, gd_index] = max(pks); % not used
%             woof = diff(tempY);
%             gd_index = find(abs(woof) == max(abs(woof)));
%             if length(gd_index) > 1
%                 gd_index = gd_index(1);
%             end
%             if gd_index > length(locs)
%                 warning('Need to figure out what happens here')
%             else
%                 gd_tp = locs(gd_index) + ss_tp + 1 - round(3/TR) - floor(w(gd_index)*1.5); %contrast arrival
%                 pk_tp = locs(gd_index) + ss_tp - round(3/TR) + 2; %time to peak
%                 
%                 
%                 h(3) = figure('Color', [1 1 1],'Visible','off');
%                 title('Mean Time Course Across Volume'); hold on;
%                 plot(tempY, 'k', 'LineWidth', 2'); hold on;
%                 if ne == 2 %LCB 4/6/2016 % ne is the number of echos NS
%                     I2 = squeeze(vol(:,:,:,1,:)); %only need first echo dataset for this section
%                     mean_tc2 = mean(reshape(I2, [nx*ny*nz nt])); %mean time course
%                     plot(mean_tc2, '--k', 'LineWidth', 2'); hold on;
%                     legend('First Echo', 'Second Echo');
%                     %         elseif ne == 1 % NS 20200402, might not work
%                     %             I2 = squeeze(vol(:,:,:,1,:)); %only need first echo dataset for this section
%                     %             mean_tc2 = mean(reshape(I2, [nx*ny*nz nt])); %mean time course
%                     %             plot(mean_tc2, '--k', 'LineWidth', 2'); hold on;
%                     %             legend('First Echo');
%                 end
%                 
%                 xlabel('Time points', 'FontSize', 14); ylabel('Signal A.U.', 'FontSize', 14);
%                 ylim = get(gca, 'ylim');
%                 line([ss_tp ss_tp], ylim, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
%                 %         line([gd_tp gd_tp], ylim, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
%                 line([gd_index gd_index], ylim, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
%                 str1 = '\rightarrow Steady State Reached';
%                 text(ss_tp, tempY(ss_tp)*1.25, str1, 'FontSize',14); %LCB 4/6/2016
%                 str1 = '\rightarrow Contrast Arrival';
%                 text(gd_tp, tempY(gd_tp)*1.1, str1,'FontSize',14); %LCB 4/6/2016
%                 %
%                 %             [envHigh, envLow] = envelope(tempY(:,ss_tp:end),16,'peak');
%                 %             envMean = (envHigh+envLow)/2;
%                 %
%                 %             plot(1:150,tempY(:,ss_tp:end), ...
%                 %                 1:150,envHigh, ...
%                 %                 1:150,envMean, ...
%                 %                 1:150,envLow)
%                 %             axis tight
%                 %             legend('Raw Signa','High','Mean','Low','location','best')
%                 %             ylabel('Signal AU')
%                 %             xlabel('Time (ms)')
%                 %             title('Signal vs Time for Voxel')
%             end
%         end
        return;
end

end