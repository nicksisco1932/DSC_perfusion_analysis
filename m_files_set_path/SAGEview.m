function SAGEview(data,varargin)
% data(1,1,:)=sin([pi/20:pi/20:20*pi]);
% data=repmat(data,[20,20,1]);
% data=data+rand(size(data));

ts=[];
%=======read in optional inputs=========
for ii=1:length(varargin)
    if isequal(varargin{ii},'ts')
        tmp=varargin{ii+1};
        tmp=(tmp-min(tmp(:)));
        tmp=tmp./max(tmp(:));
        tmp=tmp.*range(data(:));
        tmp=tmp+min(data(:));
        ts(:,ii)=cat(1,ts,tmp(:));
    elseif isequal(varargin{ii},'image')
        IM=varargin{ii+1};
    elseif isequal(varargin{ii},'colorchoice')
        colorchoice=varargin{ii+1};
    elseif isequal(varargin{ii},'TR')
        TR=varargin{ii+1};
    end
end
%=======================================
if ~exist('colorchoice','var');
    colorchoice='gray';
end
% % if ~exist('TR','var');
% %     TR=1;
% %     warning('assuming temporally sampling at 1Hz');
% % end

%transform the data to have two spatial dimensions and one temporal
%dimension
% switch ndims(data)
%     case 3
%     case 4
%         disp('montaging across the 3rd dimension of ''data''');
%         
%         dataa=layoutfcn(squeeze(data(:,:,:,1)));
%         dataa=zeros([size(dataa) size(data,4)]);
%         
%         for ii=1:size(data,4);
%             dataa(:,:,ii)=layoutfcn(squeeze(data(:,:,:,ii)));
%         end;
%             
%     otherwise
%         error('the data matrix must be 3 or 4 dimensional');
% end
% data=dataa;clear dataa;




%create an image to show if none were passed
if ~exist('IM','var');
    IM=data;
    size(IM);
end

%make sure that the IM is 2D
switch ndims(IM)
    case 2
    case 3
        disp('montaging across the 3rd dimension of ''image''');        
        IM=layoutfcn(IM);  
    case 4
        disp('Plot all 5 echoes');
        str = sprintf('\nBlue: TE1 \nGreen: TE2 \nRed: TE3\nTeal: TE4\nPurple: TE5\n');
        disp(str);
        IM = mean(abs(IM(:,:,1,:)),4);
    otherwise
        error('the IM matrix must be 2 or 3 dimensional');
end


h=figure;
set(h,'WindowButtonDownFcn',@wbdfcn);   %setup the action for the 1st button press

%----- display the mean image in the left hand plane -----
g=subplot(1,2,1);
imagesc(IM); colormap('gray')
% if ndims(IM)==2;
%     mycolormap(colorchoice);
% end;
%axis square;
axis image;
%-----

%----- open the 2nd axes -----
k=subplot(1,2,2);
axis square;
grid minor;
%-----

    function wbdfcn(src,evnt)  %by putting the 'end' at the end of each function, matlab can tell which functions are nested, and will give then access to their parent function's workspace (apparently).
        %-----make sure that only the image is contained in the left hand axes-----
        gchildren=get(g,'children');       
        for ii=1:length(gchildren);
            if ~isequal(get(gchildren(ii),'Type'),'image');
                delete(gchildren(ii));
                drawnow;
            end
        end;
        clear ii;
        %-----
        
        %----- set pointer style and get current position -----
        set(h,'Pointer','crosshair');
        tmp=round(get(g,'CurrentPoint'));
        y=tmp(1,1);
        x=tmp(1,2);
        
        %-----
        
        %----- ensure that the current position is within the current axes and if so then plot the time course for that voxel -----
        if x>0 & x<=size(data,1) & y>0 & y<=size(data,2)
            axes(k);
                ts2=(ts.*range(squeeze(data(x,y,:))))+min(squeeze(data(x,y,:)));
% %                 plot(1:size(data,4),cat(2,ts2,squeeze(data(x,y,1,:)))',1:size(data,4),cat(2,ts2,squeeze(data(x,y,2,:)))',1:size(data,4),cat(2,ts2,squeeze(data(x,y,3,:)))',1:size(data,4),cat(2,ts2,squeeze(data(x,y,4,:)))',1:size(data,4),cat(2,ts2,squeeze(data(x,y,5,:)))','Linewidth',1.5);axis square;
                plot(1:size(data,4),cat(2,ts2,squeeze(data(x,y,:,:)))','Linewidth',1.5);axis square;
                title(['your current position is: ',num2str([tmp(1,1),tmp(1,2)])]);
% %             legend('dual echo','te1','te2','BW','TRATE','Location','SouthWest');
            drawnow;
        end
        set(h,'WindowButtonMotionFcn',@wbmfcn)
        set(h,'WindowButtonDownFcn',@wbdfcn2);
        %-----
        
        function wbmfcn(src,evnt)
            %----- repeat the plotting for each new position while the pointer is in motion -----
            tmp=round(get(g,'CurrentPoint'));
            y=tmp(1,1);
            x=tmp(1,2);

            if x>0 & x<=size(data,1) & y>0 & y<=size(data,2)
                axes(k);
                try
                    ts2=(ts.*range(squeeze(data(x,y,:))))+min(squeeze(data(x,y,:)));
                    plot(1:size(data,4),cat(2,ts2,squeeze(data(x,y,:,:)))','Linewidth',1.5);axis square;
% %                     plot(1:size(data,4),cat(2,ts2,squeeze(data(x,y,1,:)))',1:size(data,4),cat(2,ts2,squeeze(data(x,y,2,:)))',1:size(data,4),cat(2,ts2,squeeze(data(x,y,3,:)))',1:size(data,4),cat(2,ts2,squeeze(data(x,y,4,:)))',1:size(data,4),cat(2,ts2,squeeze(data(x,y,5,:)))','Linewidth',1.5);axis square;
                    %                 catch
                    %                     tmp=lasterror;
                    %                     disp(tmp.message);
                end
                title(['your current position is: ',num2str([tmp(1,1),tmp(1,2)])]);%%legend('TE1','TE2','TE3','TE4','TE5');
% %                 legend('dual echo','te1','te2','BW','TRATE','Location','SouthWest');
                drawnow;
            end
            %-----
        end
        function wbdfcn2(src,evnt)
            %----- plot the current point and return to the initial state upon the 2nd button press -----
            set(h,'WindowButtonMotionFcn','');
            set(h,'WindowButtonDownFcn',@wbdfcn,'Pointer','arrow');
            
            axes(k)
            %fftanalysis(squeeze(data(x,y,:)),1./TR,[],'d');
            
            axes(g);
            hold all;
            q=plot(xlim,[x x],'g-',[y y],ylim,'g-');
            set(q,'markersize',1,'linewidth',1);
            clear q;
            hold off;
            drawnow;
            %-----
        end

    end
end


function [laidout]=layoutfcn(threed)
    sizes=size(threed);
    if length(sizes)==2
        sizes=[sizes 1];
    end
    side=ceil(sqrt(sizes(3)));
    laidout=zeros(sizes(1)*ceil(sizes(3)/side),sizes(2)*side);
    for m=1:sizes(3)
        q=m/side;
        if q==floor(q)
            i=side; %no zero padding needed
        else
            i=round((q-floor(q))*side); %zero padding needed
        end
        j=ceil(q);
        laidout([1+(sizes(1)*(j-1)):sizes(1)+(sizes(1)*(j-1))],[1+(sizes(2)*(i-1)):sizes(2)+(sizes(2)*(i-1))])=squeeze(threed(:,:,m));
    end;
end
