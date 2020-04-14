function [varargout] = custommontage(varargin)

switch nargin
    case 0
        error('No inputs passed');
    case 1
        IM_in=varargin{1};
        col=ceil(sqrt(size(IM_in(:,:,:),3)));
        row=ceil(sqrt(size(IM_in(:,:,:),3)));
    case 2
        IM_in=varargin{1};
        col=varargin{2};
        row=ceil(size(IM_in,3)./col);
    case 3
        IM_in=varargin{1};
        col=varargin{2};
        row=varargin{3};
        if (col*row)<size(IM_in(:,:,:),3)
            error('The specified matrix size is not large enough to hold all the slices contained in the input image.');
        end
    otherwise
        error('unknown inputs');
end

switch ndims(IM_in)
    case 2
        varargout{1}=IM_in;
    case 3
        sizes=size(IM_in);
        
        IM_out=[];
        counter=0;
        
        for jj=1:row
            tmp=[];
            for kk=1:col
                counter=counter+1;
                if counter<=size(IM_in,3)
                    tmp=cat(2,tmp,IM_in(:,:,counter));
                else
                    tmp=cat(2,tmp,zeros(size(IM_in,1),size(IM_in,2)));
                end
            end
            IM_out=cat(1,IM_out,tmp);
        end
        
        varargout{1}=IM_out;
    case 4
        sizes=size(IM_in);
        
        IM_out=[];
        counter=0;
        
        for jj=1:row
            tmp=[];
            for kk=1:col
                counter=counter+1;
                if counter<=size(IM_in(:,:,:),3)
                    tmp=cat(2,tmp,IM_in(:,:,counter));
                else
                    tmp=cat(2,tmp,zeros(size(IM_in,1),size(IM_in,2)));
                end
            end
            IM_out=cat(1,IM_out,tmp);
        end
        
        varargout{1}=IM_out;
        
    case 5
        sizes=size(IM_in);
        
        IM_out=[];
        counter=0;
        
        for jj=1:row
            tmp=[];
            for kk=1:col
                counter=counter+1;
                if counter<=size(IM_in(:,:,:),3)
                    tmp=cat(2,tmp,IM_in(:,:,counter));
                else
                    tmp=cat(2,tmp,zeros(size(IM_in,1),size(IM_in,2)));
                end
            end
            IM_out=cat(1,IM_out,tmp);
        end
        
        varargout{1}=IM_out;
        
    otherwise
        error('INPUT IMAGE HAS TOO MANY DIMENSIONS');
end

