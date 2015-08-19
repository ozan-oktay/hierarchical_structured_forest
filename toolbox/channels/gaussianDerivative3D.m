function gaussianDerivative3D (inputArr, outputArr)

imdata    = imread('ngc6543a.jpg');
inputArr  = repmat(imdata(:,:,1),[1,1,100]);
inputSize = size(inputArr);
inputArr  = single(inputArr);

r_val    = 5;
sigma    = 2
num_comp = numelements(r_val) * 6;
outputArr= zeros([inputSize,num_comp],'single');
filters  = zeros(2*r_val+1,2*r_val+1,2*r_val+1,num_comp);

% Generate Gaussian Derivative Filters
for loopId=1:numelements(r_val)
  index = (1:6) + (loopId-1)*6;
  filters(:,:,:,index) = generateFilters(r_val(loopId), sigma(loopId));
end

% Convolve with Gaussian Derivative Filter
for loopId=1:numelements(r_val)
  outputArr(:,:,:,loopId) = convn(inputArr,filters(:,:,:,loopId),'same');
end

% Square it and aggregate it over spatial region
outputArr = outputArr.^2;


end

function out = generateFilters(r, sigma)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[x,y,z] = ndgrid(-r:r,-r:r,-r:r);
Gk = ((sqrt(2*pi)*sigma)^3) .* exp(-(x.^2+y.^2+z.^2)/(2*sigma^2));

[Gx,Gy,Gz] = gradient(Gk);   
[Gxx,Gxy,Gxz] = gradient(Gx);
[~,Gyy,Gyz]   = gradient(Gy);
[~,~,Gzz]     = gradient(Gz);

Gxx = Gxx / sum(Gxx(:).^2);
Gxy = Gxy / sum(Gxy(:).^2);
Gxz = Gxz / sum(Gxz(:).^2);
Gyy = Gyy / sum(Gyy(:).^2);
Gyz = Gyz / sum(Gyz(:).^2);
Gzz = Gzz / sum(Gzz(:).^2);
out = cat(4,Gxx,Gxy,Gxz,Gyy,Gyz,Gzz);

end

