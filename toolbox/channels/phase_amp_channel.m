%PHASE_CHANNEL:
%This function computes feature assymetry images for given 3D input
%ultrasound image. 

%It uses local phase and local magnitude properties of the signal to
%compute feature assymetry

function [ outputImage ] = phase_amp_channel( inputImage, b_normalize, b_fs_fa, bw_factor )

if (nargin < 4)
  bw_factor = 1.25;
end

if (nargin < 3)
  b_fs_fa = 'feature_asymmetry';
end

loaddependencies();
inputType  = class(inputImage);
inputImage = double(inputImage);

mask            = myMaskGenerator(inputImage); %%% Generate Mask for boundaries
wavLengths      = [4.0,6.0,8.0] * bw_factor;  % center frequency of bandpass filter LOP
[~,FA,localEng] = phaseImageFunction(inputImage,wavLengths,'LOP',b_fs_fa); %%% Generate Phase Image 

FA (mask)   = 0.0;
FA_ENG      = FA.*localEng(:,:,:,2); %%% Perform Post-processing 
if (b_normalize)
  FA_ENG = invertImage(FA_ENG);
  FA_ENG = myNormalizeImage (FA_ENG);
end
outputImage = cast(FA_ENG, inputType);

end


function loaddependencies ()
str            = '/homes/oo2113/workspace/matlab/nifti-read/';         addpath (genpath(str));
str            = '/homes/oo2113/workspace/matlab/displayFunctions/';   addpath (genpath(str));
end

function outputImage = myNormalizeImage (inputImage)
max_val     = max(inputImage(:));
min_val     = min(inputImage(:));
slope       = 1024 / (max_val-min_val);
intersec    = -1 * slope * min_val;
outputImage = inputImage * slope + intersec;
end

function outputMask = myMaskGenerator (inputImage)
mask = (inputImage <= 0);
mask(1,:,:) = true; mask(end,:,:) = true;
mask(:,1,:) = true; mask(:,end,:) = true;
mask(:,:,1) = true; mask(:,:,end) = true;

[x,y,z] = ndgrid(-3:3); iter=2;
morpElement = strel(sqrt(x.^2 + y.^2 + z.^2) <=3);
for ii=1:iter,   mask = imerode(mask,morpElement); end;
for ii=1:iter+2, mask = imdilate(mask,morpElement); end;
outputMask   = mask;
end

function outputImage = invertImage (inputImage)
max_val     = max(inputImage(:));
min_val     = min(inputImage(:));
slope       = 1 / (max_val-min_val);
intersec    = -1 * slope * min_val;
outputImage = inputImage * slope + intersec;
outputImage = 1 - outputImage;
end

function [phaseAng,FA,localEng] = phaseImageFunction(image,wavLengths,bP_flag,b_fs_fa)

% ========================================================================
% Copyright(c) 
% All Rights Reserved.
% ----------------------------------------------------------------------

%In fact, the Riesz kernels can also be pre-computed.
[rows, cols, slices] = size (image);

if (strcmp(bP_flag,'LOP'))
    [Rx, Ry, Rz, BP] = generateRieszKernel_3D_LOP(rows, cols, slices, wavLengths);
elseif (strcmp(bP_flag,'GD'))
    [Rx, Ry, Rz, BP] = generateRieszKernel_3D_GD (rows, cols, slices, wavLengths);
else
    error ('Undefined Bandpass Filter');
end

%the image need some preprocessing
image    = double(image);
Ex       = mean2(image);
sigma    = std(image(:));
image    = (image - Ex) ./ sigma;
FA       = zeros (size(image));
phaseAng = zeros (cat(2,size(image),numel(wavLengths)));
localEng = zeros (cat(2,size(image),numel(wavLengths)));

for scaleIndex = 1:numel(wavLengths)
    
    fftIm = fftn(double(image));
    bpIm  = real(ifftn(fftIm .* BP(:,:,:,scaleIndex))); 
    RxRes = real(ifftn(fftIm .* Rx(:,:,:,scaleIndex)));
    RyRes = real(ifftn(fftIm .* Ry(:,:,:,scaleIndex)));
    RzRes = real(ifftn(fftIm .* Rz(:,:,:,scaleIndex)));
      
    localEng(:,:,:,scaleIndex) = sqrt(RxRes.^2 + RyRes.^2 + RzRes.^2 + bpIm.^2);
    phaseAng(:,:,:,scaleIndex) = atan2(sqrt(RxRes .^ 2 + RyRes .^ 2 + RzRes .^ 2), bpIm);
    
    even  = bpIm;
    odd   = sqrt(RxRes .^ 2 + RyRes .^ 2 + RzRes .^ 2);
    T     = exp(mean2(log(sqrt(even.^2 + odd.^2))));
    if (strcmp(b_fs_fa,'feature_asymmetry')), upper = (abs(odd) - abs(even)) - T;
    elseif (strcmp(b_fs_fa,'feature_symmetry')), upper = (abs(even) - abs(odd)) - T;end
    FA    = FA + ( (upper.*heaviside(upper)) ./ (sqrt(even.^2 + odd.^2) + 1e-9) ); 

end
end

%==========================================================================
%This function returns the Riesz transforms kernels and the band-pass
%filter kernel (LOP)

function [Rx, Ry, Rz, LOP] = generateRieszKernel_3D_LOP (rows, cols, slices, wavLengths)

    [u1, u2, u3] = meshgrid( ([1:cols]-(fix(cols/2)+1))/(cols-mod(cols,2)), ...
                         ([1:rows]-(fix(rows/2)+1))/(rows-mod(rows,2)), ...
                         ([1:slices]-(fix(slices/2)+1))/(slices-mod(slices,2)) );

    u1 = ifftshift(u1);  
    u2 = ifftshift(u2);
    u3 = ifftshift(u3);
    
    radius        = sqrt(u1.^2 + u2.^2 + u3.^2);    
    radius(1,1,1) = 1;
    
    R1 = -1i*u1./radius;  
    R2 = -1i*u2./radius;
    R3 = -1i*u3./radius;
    radius(1,1,1) = 0;
    
    for wavIndex = 1:length(wavLengths)
        LOP(:,:,:,wavIndex) = -8  * pi^3 * (radius .^2) .* exp(-2*pi * radius * wavLengths(wavIndex));
         Rx(:,:,:,wavIndex) = R1 .* LOP(:,:,:,wavIndex); 
         Ry(:,:,:,wavIndex) = R2 .* LOP(:,:,:,wavIndex); 
         Rz(:,:,:,wavIndex) = R3 .* LOP(:,:,:,wavIndex);
    end
end

function [Rx, Ry, Rz, GD] = generateRieszKernel_3D_GD (rows, cols, slices, sigma)


    [u1, u2, u3] = meshgrid( ([1:cols]-(fix(cols/2)+1))/(cols-mod(cols,2)), ...
                         ([1:rows]-(fix(rows/2)+1))/(rows-mod(rows,2)), ...
                         ([1:slices]-(fix(slices/2)+1))/(slices-mod(slices,2)) );

    u1 = ifftshift(u1);  
    u2 = ifftshift(u2);
    u3 = ifftshift(u3);
    
    radius        = sqrt(u1.^2 + u2.^2 + u3.^2);    
    radius(1,1,1) = 1;
    
    R1 = -1i*u1./radius;  
    R2 = -1i*u2./radius;
    R3 = -1i*u3./radius;
    radius(1,1,1) = 0;
    
    for sigmaInd = 1:length(sigma)
        GD(:,:,:,sigmaInd)  = radius .* exp (-1 * radius.^2 * sigma(sigmaInd)^2);
         Rx(:,:,:,sigmaInd) = R1 .* GD(:,:,:,sigmaInd); 
         Ry(:,:,:,sigmaInd) = R2 .* GD(:,:,:,sigmaInd); 
         Rz(:,:,:,sigmaInd) = R3 .* GD(:,:,:,sigmaInd);
    end
end