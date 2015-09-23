function [chnsReg,chnsSim] = edgesChns_3D(I, opts)

% Compute features for structured edge detection.
%
% For an introductory tutorial please see edgesDemo.m.
%
% USAGE
%  [chnsReg,chnsSim] = edgesChns( I, opts )
%
% INPUTS
%  I          - [h x w x d] color input image
%  opts       - structured edge model options
%
% OUTPUTS
%  chnsReg    - [h x w x d x nChannel] regular output channels
%  chnsSim    - [h x w x d x nChannel] self-similarity output channels
%
%
% Original version: 
%
% Structured Edge Detection Toolbox      Version 3.01
% Code written by Piotr Dollar, 2014.
%
% Updated version: 
% 
% Author: Ozan Oktay
% Date:   May 2015
% Email:  o.oktay13@imperial.ac.uk
%

shrink=opts.shrink; chns=cell(1,opts.nChns); k=0;
Ishrink=imResample3D(I,size(I)/shrink); k=k+1; chns{k}=Ishrink;  

for i = 1:2, s=2^(i-1);
if(s==shrink), I1=Ishrink; else I1=imResample3D(I,size(I)/s); end
I1 = convTri3D( I1, opts.grdSmooth );
M  = gradientMag3D( I1, opts.normRad, 0.05 ); 
k=k+1; chns{k}=imResample3D(M,size(M)*s/shrink);
end

chndim = ndims(I)+1;
imdim  = ndims(I);
chns   = cat(chndim,chns{1:k}); assert(size(chns,chndim)==opts.nChns);
chnSm  = opts.chnSmooth/shrink; if(chnSm>1), chnSm=round(chnSm); end
simSm  = opts.simSmooth/shrink; if(simSm>1), simSm=round(simSm); end

chnsReg=convTri3D(chns,chnSm,[ones(1,imdim),0]); 
chnsSim=convTri3D(chns,simSm,[ones(1,imdim),0]);

end