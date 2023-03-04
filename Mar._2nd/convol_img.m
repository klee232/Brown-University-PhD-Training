function conv_img=convol_img(in_img,kernel_size,varargin)
% get the dimension of the input image
wid_img=size(in_img,1);
len_img=size(in_img,2);
chn_img=size(in_img,3);

% default setting for convolution
defaultkernel=rand(kernel_size,kernel_size,chn_img);
defaultstride=1;
defaulttype="valid";
% check input functions
checkarray = @(x) ~isempty(x);
checkinteger = @(x) isfinite(x) && x == floor(x);
% insert input parser for organizing the required and optional input
p=inputParser;
addRequired(p,'in_img',checkarray);
addRequired(p,'kernel_size',checkinteger);
addParameter(p,'kernel',defaultkernel,checkarray);
addParameter(p,'stride',defaultstride,checkinteger);
addParameter(p,'type',defaulttype);
parse(p,in_img,kernel_size,varargin{:});
% store the value to assign variable
in_img=p.Results.in_img;
kernel_size=p.Results.kernel_size;
kernel=p.Results.kernel;
stride=p.Results.stride;
type=p.Results.type;
% check if the kernel is a single channel kernel
if size(kernel,3)==1
   kernel=cat(3,kernel,kernel,kernel);
end


% allocate array and pad size according convolution type
if strcmp(type,"valid")
    pad_width=0;
    pad_height=0;
    img_temp=in_img;
    wid_conv=floor(double(wid_img+2*pad_width-kernel_size)/double(stride))+1;
    hei_conv=floor(double(len_img+2*pad_height-kernel_size)/double(stride))+1;
    chn_conv=size(img_temp,3);
    conv_img=zeros(wid_conv,hei_conv,chn_img);
elseif strcmp(type,"same")
    pad_width=(kernel_size-1)/2;
    pad_height=(kernel_size-1)/2;
    pad_hor=zeros(pad_width,len_img,chn_img);
    img_temp=cat(1,pad_hor,in_img,pad_hor);
    wid_temp=size(img_temp,1);
    pad_ver=zeros(wid_temp,pad_height,chn_img);
    img_temp=cat(2,pad_ver,img_temp,pad_ver);
    wid_conv=double(wid_img+2*pad_width-kernel_size)/double(stride)+1;
    hei_conv=double(len_img+2*pad_height-kernel_size)/double(stride)+1;
    chn_conv=size(img_temp,3);
    conv_img=zeros(wid_img,len_img,chn_img);
else
    error("something with the type setting. please check");
end

% conducting convolution on the input image
for i_wid=1+floor(kernel_size/2):stride:wid_conv-floor(kernel_size/2)
    for i_col=1+floor(kernel_size/2):stride:hei_conv-floor(kernel_size/2)
        for i_chn=1:chn_conv
            % region of interest processing
            region_of_interest=img_temp(i_wid-floor(kernel_size/2):i_wid+floor(kernel_size/2),...
                                        i_col-floor(kernel_size/2):i_col+floor(kernel_size/2),...
                                        i_chn);
            % convolution part
            % conduct convolution (be sure to convert the original datatype
            % to double to conduct the calculation)
            conv_region_of_interest=double(region_of_interest).*double(kernel(:,:,i_chn));
            conv_region_of_interest=sum(conv_region_of_interest,'all');
            % store the outcome into assigned position
            conv_img(i_wid,i_col,i_chn)=conv_region_of_interest;
        end
    end
end


end