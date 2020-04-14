function dR2s = process_deltaR2star(vol, TIMEarray, TParray, varargin)
%process_deltaR2star This function determines the delta R2 star curves
%
% Laura Bell 1/5/2016 (Code from Ashley Stokes)
%
% Usage: dR2s = process_R2star(vol, TIMEarray, TParray, varargin)  
%        dR2s = process_R2star(vol, TIMEarray, TParray, 'LeakageCorrection', 'Boxerman')
%        dR2s = process_R2star(vol, TIMEarray, TParray, 'LeakageCorrection', 'Bjornerud', AIF_dR2star)
% Inputs:                                                                                                                  
% - vol is a 4D [nx ny nz nt] or 5D [nx ny nz ne nt] matrix - recommended to be masked already                            
% - TIMEarray is an array with TR and then TEs in seconds - shortest to longest TE                                         
% - TParray is an array with time points specified for baseline calculations 
% - Optional inputs:
%       1) 'Boxerman' 
%       2) 'Bjornerud' - array of R2star values for the AIF is needed for the Bjornerud correction method
%
% Additional functions needed:
% - produceEVthreshold: gets called for Bjornerud method
%                                                                                                                                                                                                               
% Outputs:                                                                                                                 
% - output will a structure of delta R2 star curve array, if multiple echos acquired then multiple R2 star curves produced            

%% Get optional inputs
if length(varargin) >= 1
    if isequal(varargin{1},'Boxerman')
        LCflag = 'Boxerman';
        display(sprintf('\nAssuming %s leakage correction method for delta R2 star curves.', LCflag));
    elseif isequal(varargin{1}, 'Bjornerud')
        LCflag = 'Bjornerud';
        AIF = varargin{2};
        display(sprintf('\nAssuming %s leakage correction method for delta R2 star curves.', LCflag));
    else
        warning('Leakage correction method not recognized.');
        return;
    end
end

if ~exist('LCflag', 'var')
    LCflag = 'Uncorrected';
    display(sprintf('\nAssuming %s method for delta R2 star curves.', LCflag));
end

%% Read in volume and determine matrix size
if length(size(vol)) == 2
    [nx, nt] = size(vol); %nx is number of pixels; nt is number of time points
    ny = 1; nz = 1;
    display(sprintf('Number of pixels: %d, Number of timepoints: %d \n', nx, nt));
    echoFlag = 0;
    arrayFlag = 1;
    
    I_v = vol;
    I_v(I_v == 0) = nan;
    TE1 = TIMEarray(2);
    TR = TIMEarray(1);
elseif length(size(vol)) == 5
    [nx,ny,nz,ne,nt] = size(vol);
    warning('This code only uses information for 2 echoes.'); %you'll need to code if you want to add other echoes
    display(sprintf('\nThis is a multi echo dataset with %d echoes: ', ne));
    display(sprintf('nx: %d ny: %d nz: %d nt: %d \n', nx, ny, nz, nt));
    echoFlag = 1;
    if ne>2
        seFlag = 1;
        TEse = TIMEarray(end);
    else
        seFlag = 0;
    end
    
    I_v = squeeze(reshape(vol, [nx*ny*nz ne nt]));
    I_v(I_v == 0) = nan;
    TE1 = TIMEarray(2);
    TE2 = TIMEarray(3);
    TR = TIMEarray(1);
elseif length(size(vol)) == 4
    [nx,ny,nz,nt] = size(vol);
    display(sprintf('\nThis is a single echo dataset: '));
    display(sprintf('nx: %d ny: %d nz: %d nt: %d \n', nx, ny, nz, nt));
    echoFlag = 0;seFlag = 0;
    
    I_v = reshape(vol, [nx*ny*nz nt]);
    I_v(I_v == 0) = nan;
    TE1 = TIMEarray(2);
    TR = TIMEarray(1);
end

if exist('arrayFlag', 'var') && ~strcmp(LCflag, 'Uncorrected')
    warning('Can not have an array for input for leakage corrected methods. Ending script.');
    return;
end
%% set time points
ss_tp = TParray(1); %steadystate reached (1 if no dummy scans included)
gd_tp = TParray(2); %gad arrival time

%% uncorrected R2star curves
switch echoFlag
    case 0
        S0_TE1 = squeeze(nanmean(I_v(:,ss_tp:gd_tp),2));
        S0_TE1 = repmat(S0_TE1, [1 nt]);
        
        dR2s_TE1 = log(S0_TE1./squeeze(I_v(:,:)))/TE1;
        dR2s_TE1(isinf(dR2s_TE1)) = 0;
        dR2s_TE1(isnan(dR2s_TE1)) = 0;
    case 1
        
            S0_TE1 = squeeze(nanmean(I_v(:,1,ss_tp:gd_tp),3)); S0_TE1 = repmat(S0_TE1, [1 nt]);
            S0_TE2 = squeeze(nanmean(I_v(:,2,ss_tp:gd_tp),3)); S0_TE2 = repmat(S0_TE2, [1 nt]);

            dR2s_TE1 = log(S0_TE1./squeeze(I_v(:,1,:)))/TE1; dR2s_TE1(isinf(dR2s_TE1)) = 0; dR2s_TE1(isnan(dR2s_TE1)) = 0;
            dR2s_TE2 = log(S0_TE2./squeeze(I_v(:,2,:)))/TE2; dR2s_TE2(isinf(dR2s_TE2)) = 0; dR2s_TE2(isnan(dR2s_TE2)) = 0;
            dR2s_all = (1/(TE2-TE1)).*(log( (squeeze(I_v(:,1,:))./S0_TE1) ./ (squeeze(I_v(:,2,:))./S0_TE2) ));
            dR2s_all(isinf(dR2s_all)) = 0; dR2s_all(isnan(dR2s_all)) = 0; 
        
        switch seFlag
            case 1
                STE0 = squeeze(I_v(:,1,:).*(I_v(:,1,:)./I_v(:,2,:)).^(TE1/(TE2-TE1)));
                S0_TE0 = squeeze(nanmean(STE0(:,ss_tp:gd_tp),2)); S0_TE0 = repmat(S0_TE0, [1 nt]);
                S0_TEse = squeeze(nanmean(I_v(:,end,ss_tp:gd_tp),3)); S0_TEse = repmat(S0_TEse, [1 nt]);
                
                dR2s_TEse = log(S0_TEse./squeeze(I_v(:,end,:)))/TEse; dR2s_TEse(isinf(dR2s_TEse)) = 0; dR2s_TEse(isnan(dR2s_TEse)) = 0;
                dR2s_allSE = (1/(TEse)).*(log( (squeeze(STE0(:,:))./S0_TE0) ./ (squeeze(I_v(:,end,:))./S0_TEse) ));
                dR2s_allSE(isinf(dR2s_allSE)) = 0; dR2s_allSE(isnan(dR2s_allSE)) = 0;
        end
end

%% Start leakage correction if requested
switch LCflag
    case 'Boxerman'
        switch echoFlag
            case 0
                mnSI_end = nanmean(I_v(:,end-10:end),2);
                pos_threshold = nanmean(I_v(:,ss_tp:gd_tp),2) + nanstd(squeeze(I_v(:,ss_tp:gd_tp)), [], 2);
                neg_threshold = nanmean(I_v(:,ss_tp:gd_tp),2) - nanstd(squeeze(I_v(:,ss_tp:gd_tp)), [], 2);
                dR2s_v = dR2s_TE1;
            case 1
                mnSI_end = nanmean(I_v(:,2,end-10:end),3);
                pos_threshold = nanmean(I_v(:,2,ss_tp:gd_tp),3) + nanstd(squeeze(I_v(:,2,ss_tp:gd_tp)), [], 2);
                neg_threshold = nanmean(I_v(:,2,ss_tp:gd_tp),3) - nanstd(squeeze(I_v(:,2,ss_tp:gd_tp)), [], 2);
                dR2s_v = dR2s_all;
        end
        
        % First calculate non enhancing voxels
        nonenhance_map1 = zeros(nx,ny,nz);
        nonenhance_map2 = zeros(nx,ny,nz);
        
        nonenhance_map1(mnSI_end < pos_threshold) = 1;
        nonenhance_map2(mnSI_end > neg_threshold) = 1;
        nonenhance_map = nonenhance_map1 + nonenhance_map2;
        nonenhance_map(nonenhance_map == 1) = 0;
        nonenhance_map(nonenhance_map == 2) = 1;
        %figure, montage(permute(nonenhance_map, [1 2 4 3]));
        %save nonenhance_map nonenhance_map
        
        % Average delta R2 star based on non-enhancing mask
        dR2s_NEvec = dR2s_v .* repmat(nonenhance_map(:),[1 nt]);
        dR2s_NEvec(dR2s_NEvec == 0) = NaN;
        for x = 1:length(dR2s_NEvec)
            if ~isreal(dR2s_NEvec(x,:))
                dR2s_NEvec(x,:) = nan;
            end
        end
        dR2s_WBNE = nanmean(dR2s_NEvec);
        baseline_dR2sWBNE = nanmean(dR2s_WBNE(ss_tp:gd_tp));
        dR2s_WBNE(1:ss_tp) = repmat(baseline_dR2sWBNE, [1 ss_tp]);
        
        %         h = figure;
        %         h.Visible = 'on';
        %         plot(dR2s_WBNE, 'LineWidth', 2);
        %         title('Whole brain non-enhancing ROI');
        %         xlabel('Time points'); ylabel('delta R2 star');
        %         savefig(h, 'Boxerman_WholeBrainNonEnhancing.fig');
        
        % Determine K1, K2, and corrected delta R2star curves
        K1 = zeros(nx*ny*nz,1); %scales blood volume
        K2 = zeros(nx*ny*nz,1); %scales permeability
        dR2s_Boxerman = zeros(nx*ny*nz,nt);
        
        x0 = [1 0];
        options = optimset('TolFun',1e-9,'Tolx',1e-9,'MaxIter',1000,'Display','off');
        int_dR2sWBNE = cumsum(dR2s_WBNE);
        
        display(sprintf('\n About to start Boxerman K1 and K2 fitting . . . \n'));
        progressbar = waitbar(0,'Boxerman stuff going on');
        for x = 1:nx*ny*nz
            if dR2s_v(x,1) ~= 0 %assuming volume has been masked prior
                %display(['Running Boxerman fitting: ' num2str((x/(nx*ny*nz))*100) '% done.']);
                tmp_dR2s = dR2s_v(x,ss_tp:end);
                lcorrpar = lsqcurvefit('BW_lcorr',x0,1:nt-ss_tp+1,tmp_dR2s,[],[],options,dR2s_WBNE(ss_tp:end),int_dR2sWBNE(ss_tp:end));
                K1(x,1) = lcorrpar(1);
                K2(x,1) = lcorrpar(2);
                leakage = -K2(x,1).*int_dR2sWBNE(ss_tp:end);
                dR2s_Boxerman(x,ss_tp:end) = tmp_dR2s - leakage;
%                 plot(tmp_dR2s-leakage, 'r', 'LineWidth', 2), hold on, plot(tmp_dR2s); hold off;
%                 pause;
            end
            waitbar(x/(nx*ny*nz));
        end
        waitbar(1);
        close(progressbar)
        
        switch seFlag
            case 1
                mnSI_end = nanmean(I_v(:,end,end-10:end),3);
                pos_threshold = nanmean(I_v(:,end,ss_tp:gd_tp),3) + nanstd(squeeze(I_v(:,end,ss_tp:gd_tp)), [], 2);
                neg_threshold = nanmean(I_v(:,end,ss_tp:gd_tp),3) - nanstd(squeeze(I_v(:,end,ss_tp:gd_tp)), [], 2);
                dR2s_v = dR2s_allSE;
        
                % First calculate non enhancing voxels
                nonenhance_map1 = zeros(nx,ny,nz);
                nonenhance_map2 = zeros(nx,ny,nz);

                nonenhance_map1(mnSI_end < pos_threshold) = 1;
                nonenhance_map2(mnSI_end > neg_threshold) = 1;
                nonenhance_mapSE = nonenhance_map1 + nonenhance_map2;
                nonenhance_mapSE(nonenhance_mapSE == 1) = 0;
                nonenhance_mapSE(nonenhance_mapSE == 2) = 1;
 
                % Average delta R2 star based on non-enhancing mask
                dR2s_NEvec = dR2s_v .* repmat(nonenhance_mapSE(:),[1 nt]);
                dR2s_NEvec(dR2s_NEvec == 0) = NaN;
                for x = 1:length(dR2s_NEvec)
                    if ~isreal(dR2s_NEvec(x,:))
                        dR2s_NEvec(x,:) = nan;
                    end
                end
                dR2s_WBNEse = nanmean(dR2s_NEvec);
                baseline_dR2sWBNE = nanmean(dR2s_WBNEse(ss_tp:gd_tp));
                dR2s_WBNEse(1:ss_tp) = repmat(baseline_dR2sWBNE, [1 ss_tp]);

                % Determine K1, K2, and corrected delta R2star curves
                K1SE = zeros(nx*ny*nz,1); %scales blood volume
                K2SE = zeros(nx*ny*nz,1); %scales permeability
                dR2s_BoxermanSE = zeros(nx*ny*nz,nt);

                x0 = [1 0];
                options = optimset('TolFun',1e-9,'Tolx',1e-9,'MaxIter',1000,'Display','off');
                int_dR2sWBNEse = cumsum(dR2s_WBNEse);

                display(sprintf('\n About to start Boxerman K1 and K2 fitting . . . \n'));
                progressbar(0);
                for x = 1:nx*ny*nz
                    if dR2s_v(x,1) ~= 0 %assuming volume has been masked prior
                        %display(['Running Boxerman fitting: ' num2str((x/(nx*ny*nz))*100) '% done.']);
                        tmp_dR2s = dR2s_v(x,ss_tp:end);
                        lcorrpar = lsqcurvefit('BW_lcorr',x0,1:nt-ss_tp+1,tmp_dR2s,[],[],options,dR2s_WBNEse(ss_tp:end),int_dR2sWBNEse(ss_tp:end));
                        K1SE(x,1) = lcorrpar(1);
                        K2SE(x,1) = lcorrpar(2);
                        leakage = -K2SE(x,1).*int_dR2sWBNEse(ss_tp:end);
                        dR2s_BoxermanSE(x,ss_tp:end) = tmp_dR2s - leakage;
                    end
                    progressbar(x/(nx*ny*nz));
                end
                progressbar(1);
                close(progressbar)
        end
        
    case 'Bjornerud'
        switch echoFlag
            case 0
                dR2s_v = dR2s_TE1;
            case 1
                dR2s_v = dR2s_all;
        end
        
        % remove time points where steady state hasn't been reached yet and replace to keep length the same
        for ii = 1:length(dR2s_v)
            tmpdR2s = dR2s_v(ii,:);
            fill = repmat(dR2s_v(ii,ss_tp), 1, ss_tp);
            tmpdR2s(1, 1:ss_tp) = fill;
            dR2s_v(ii,:) = tmpdR2s;
        end
        
        %start Eigenvalue thresholding
        threshold = produceEVthreshold(vol, [ss_tp gd_tp]);
        
        % Do SVD - make sure masked if what to speed up process
        PerfusionMetrics = producePerfusionMetrics(AIF, dR2s_v, TR, threshold, 'AIFmatrix', 'circ');
        CBF = PerfusionMetrics.CBF;
        R_t = PerfusionMetrics.residueFunction;
        
        % Start Lorentzian fitting
        H_BL = reshape(R_t, [nx*ny*nz nt]);
        
        options = optimset('TolFun',1e-12,'LargeScale','on','Display','off','MaxIter',400);
        x0 = 2;
        
        Tcmap = zeros(nx*ny*nz,1); %Tc = capillary transit time
        [~,ind] = max(H_BL(:,1:gd_tp),[],2);
        imax = ind+4;
        
        display(sprintf('\n About to start fitting for Bjorenerud method . . . \n'));
        progressbar(0);
        for ii = 1:length(H_BL)
            tmpHBL = H_BL(ii,1);
            if tmpHBL > 0
                t = ind(ii,1):imax(ii,1);
                H_BLn=squeeze(H_BL(ii,t(1):t(end)))./max(H_BL(ii,t(1):t(end)));
                Tcpar=lsqcurvefit('lorfit',x0,t,H_BLn,0,[],options,H_BLn);
                Tcmap(ii)=Tcpar;
            end
            progressbar(ii/length(H_BL));
        end
        Tcmap(Tcmap>100) = 0;
        progressbar(1);
        close(progressbar)
        
        ka_nanmeanK = zeros(nx*ny*nz,1); %Ka is the apparent transfer constant, determined from tail end
        index = zeros(nx*ny*nz,1);
        for ii = 1:length(ka_nanmeanK)
            if Tcmap(ii,1) > 0
                index(ii,1) = ceil(Tcmap(ii,1));
                %index(ii,1) = ceil(Tcmap(ii,1)*1.5); %Ashley's code has 1.5 but I think it's TC = 1.5 times MTT (Skinner)
                ka_nanmeanK(ii,1) = nanmean(H_BL(ii,index(ii,1):end)); %1/seconds
            end
        end
        
        exC = ka_nanmeanK.*(TR).*(nt-ceil(Tcmap));
        exC(isnan(exC)) = 0;
        exC = reshape(exC, [nx ny nz]);
        
end

%% outputs
if exist('dR2s_TE1', 'var')
    dR2s.TE1 = squeeze(reshape(dR2s_TE1, [nx ny nz nt]));
end

if exist('dR2s_TE2', 'var')
    dR2s.TE2 = squeeze(reshape(dR2s_TE2, [nx ny nz nt]));
    dR2s.all = squeeze(reshape(dR2s_all, [nx ny nz nt]));
end

if exist('dR2s_TEse', 'var')
    dR2s.TEse = squeeze(reshape(dR2s_TEse, [nx ny nz nt]));
    dR2s.allSE = squeeze(reshape(dR2s_allSE, [nx ny nz nt]));
end

switch LCflag
    case 'Boxerman'
        dR2s.Boxerman = squeeze(reshape(dR2s_Boxerman, [nx ny nz nt]));
        dR2s.K1 = squeeze(reshape(K1, [nx ny nz]));
        dR2s.K2 = squeeze(reshape(K2, [nx ny nz]));
        dR2s.WBNE = squeeze(dR2s_WBNE);
        dR2s.sumdR2sWBNE = squeeze(int_dR2sWBNE);
        dR2s.nonenhancemask = squeeze(nonenhance_map);
        if seFlag
            dR2s.BoxermanSE = squeeze(reshape(dR2s_BoxermanSE, [nx ny nz nt]));
            dR2s.K1SE = squeeze(reshape(K1SE, [nx ny nz]));
            dR2s.K2SE = squeeze(reshape(K2SE, [nx ny nz]));
            dR2s.WBNEse = squeeze(dR2s_WBNEse);
            dR2s.sumdR2sWBNEse = squeeze(int_dR2sWBNEse);
            dR2s.nonenhancemaskSE = squeeze(nonenhance_mapSE);
        end
    case 'Bjornerud'
        dR2s.exC = exC;
        dR2s.CBF = CBF;
        dR2s.residueFunction = R_t;
end

end

