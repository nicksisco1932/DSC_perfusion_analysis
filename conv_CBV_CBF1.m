function [Data,Parms,AIF,AIFMask,dR2s,Perfusion] = conv_CBV_CBF1(dsc_data_orig)
Data = dsc_data_orig;

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
