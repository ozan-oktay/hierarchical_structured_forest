function temporalPCA( inputimgname, outputimgname )
%TEMPORALPCA Summary of this function goes here

% Add the following directories;
addpath(genpath('/homes/oo2113/workspace/matlab/manifold/utils'));
addpath(genpath('/vol/medic02/users/oo2113/str_forest/toolbox'));

% Read input 2DT image
inputimg   = load_untouch_nii(inputimgname);
imageSize  = size(inputimg.img);
outputType = 'single';

% Parameters
s_patch = 3;
t_patch = 2 * round(imageSize(4)/4) + 1;
n_compo = 3;
tstart  = 1;

% Output image
output = zeros(imageSize(1),imageSize(2),imageSize(3),n_compo,'single');

% Pick the patches
[X,mask] = volt2col(inputimg.img,[s_patch,s_patch,1,t_patch],tstart);
X = single(X);

% Perform pca on them - cols of U are normalized
[U,mu,~] = pca(X);
Y = pcaApply( X, U, mu, n_compo+1 );

for modeId = 1:n_compo
    output (:,:,:,modeId) = col2vol(Y(modeId+1,:), mask, imageSize(1:3),[1,1,1]);
end
  
% Rescale the output
%max_val  = max(output(:));
%min_val  = min(output(:));
%slope    = 1024 / (max_val-min_val);
%intersec = -1 * slope * min_val;
%output   = output * slope + intersec;

% Normalize the power of the transformed signal
temp_output = (output.*output);
output      = 100*(output ./ sqrt(sum(temp_output(:))));

% Save the result
inputimg.hdr.dime.dim      = [4,imageSize(1),imageSize(2),imageSize(3),n_compo,0,0,0];
inputimg.hdr.dime.datatype = 16;
inputimg.img               = cast(output,outputType);
save_untouch_nii(inputimg, outputimgname);

end

