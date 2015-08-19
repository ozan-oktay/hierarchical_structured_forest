function [E,inds,HV,pVal] = edgesDetect_3D( I, model, varargin )
% Detect edges in image.
%
% USAGE
%  [E,inds] = edgesDetect_3D( I, model )
%
% INPUTS
%  I          - [h x w x d] three dimensional input image
%  model      - structured edge model trained with edgesTrain
%
% OUTPUTS
%  E          - [h x w x d] edge probability map
%  inds       - [h/s x w/s x d/s x nTreesEval] leaf node indices
%  HV         - [h x w x d x l] hough votes for landmarks
%  pVal       - [p x 2] estimated pose values
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

% get parameters
opts=model.opts; opts.nTreesEval=min(opts.nTreesEval,opts.nTrees);
if(~isfield(opts,'sharpen')), opts.sharpen=0; end
if(~isfield(model,'segs')), model.segs=[]; model.nSegs=[]; end
opts.stride=max(opts.stride,opts.shrink); model.opts=opts;

% get other features 
stageId = model.opts.stageId;
if ( stageId == 1)
  if (nargin <= 2), error ('Expecting shape features - pem and vtk\n'); end
  oFtrs           = varargin{1};
  opts.multiscale = 0; 
end
  
if( opts.multiscale )
  % if multiscale run edgesDetect multiple times (x1 x2 resolution levels)
  ss=2.^(-1:0); k=length(ss); inds=cell(1,k); 
  siz=size(I); model.opts.multiscale=0; model.opts.nms=0; E=0;
  for i=1:k, s=ss(i); I1=imResample3D(I,siz*s);                            
    if (nargout<3 || s~=1), [E1,inds{i}]=edgesDetect_3D(I1,model);                               
    else [E1,inds{i}, HV, pVal]=edgesDetect_3D(I1,model); end
    E=E+imResample3D(E1,siz);                                                             
  end; E=E/k; model.opts=opts;
  
else
  % pad image, making divisible by 4
  siz=size(I); r=opts.imWidth/2; p=[r r r];                                
  I = padarray(I,p,'symmetric','pre');                                     
  p = p+mod(4-mod(siz(1:3)+2*r,4),4);                                      
  I = padarray(I,p,'symmetric','post');                                    
   
  % compute features and apply forest to image
  [chnsReg,chnsSim] = edgesChns_3D( I,opts );  
  if (stageId == 1), chnsShp = edgesShp_3D(oFtrs.pem,oFtrs.vtk,oFtrs.w2i,opts); else chnsShp=[]; end;
  s=opts.sharpen; if(s), I=convTri3D(I,1); end    
  
  if (nargout >= 3)
    [E,inds,HV,pVal] = houghVotes_3D_Mex(model,I,chnsReg,chnsSim,chnsShp);
    r=opts.gtWidth/2; 
    HV=HV(1+r:siz(1)+r,1+r:siz(2)+r,1+r:siz(3)+r,:); 
    HV=HV/opts.nTreesEval;
    HV=convTri3D(HV,2);
  else
    [E,inds] = houghVotes_3D_Mex(model,I,chnsReg,chnsSim,chnsShp);
  end;
  
  % normalize and finalize edge maps (each voxel eval: w^3*T/s^3)
  t=opts.stride^3/opts.gtWidth^3/opts.nTreesEval; r=opts.gtWidth/2;        
  if(s==0), t=t*2; elseif(s==1), t=t*1.8; else t=t*1.66; end
  E=E(1+r:siz(1)+r,1+r:siz(2)+r,1+r:siz(3)+r,:)*t; E=convTri3D(E,1);  
  
end
end
