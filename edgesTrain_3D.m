function model = edgesTrain_3D( varargin )
% Train structured edge detector.
%
% For an introductory tutorial please see edgesTestingDemo.m.
%
% USAGE
%  opts = edgesTrain()
%  model = edgesTrain( opts )
%
% INPUTS
%  opts       - parameters (struct or name/value pairs)
%   (1) model parameters:
%   .imWidth    - [24]  width of image patches                               
%   .gtWidth    - [12]  width of ground truth patches                        
%   (2) tree parameters:
%   .nPos       - [7e5] number of positive patches per tree
%   .nNeg       - [7e5] number of negative patches per tree
%   .nImgs      - [40]  maximum number of images to use for training
%   .nTrees     - [16]  number of trees in forest to train
%   .fracFtrs   - [1/8] fraction of features to use to train each tree      
%   .minCount   - [1]   minimum number of data points to allow split
%   .minChild   - [8]   minimum number of data points allowed at child nodes
%   .maxDepth   - [64]  maximum depth of tree
%   .discretize - ['kmeans'] options include 'pca' and 'kmeans'
%   .nSamples   - [600] number of samples for clustering structured labels 
%   .nClasses   - [2]   number of classes (clusters) for binary splits
%   .split      - ['entropy'] options include 'gini', 'entropy' and 'twoing'
%   (3) feature parameters:
%   .nOrients   - [0]   number of orientations per gradient scale             
%   .grdSmooth  - [2]   radius for image gradient smoothing (using convTri)   
%   .chnSmooth  - [2]   radius for reg channel smoothing (using convTri)
%   .simSmooth  - [4]   radius for sim channel smoothing (using convTri)
%   .normRad    - [4]   gradient normalization radius (see gradientMag)
%   .shrink     - [2]   amount to shrink channels
%   .nCells     - [4]   number of self similarity cells                       
%   (4) detection parameters (can be altered after training):
%   .stride     - [2]   stride at which to compute edges
%   .multiscale - [0]   if true run multiscale edge detector
%   .sharpen    - [2]   sharpening amount (can only decrease after training)
%   .nTreesEval - [10]  number of trees to evaluate per location
%   .nThreads   - [8]   number of threads for evaluation of trees
%   .nms        - [0]   if true apply non-maximum suppression to edges
%   (5) other parameters:
%   .seed       - [1]   seed for random stream (for reproducibility)
%   .useParfor  - [0]   if true train trees in parallel (memory intensive)
%   .modelDir   - ['models/'] target directory for storing models
%   .modelFnm   - ['ctmodel'] model filename
%   .imageDir   - ['ct_training_data/'] location of image dataset     
%   .ctmaxval   - [5000] maximum allowed intensity value - for linear scaling. 
%
% OUTPUTS
%  model      - trained structured edge detector w the following fields
%   .opts       - input parameters and constants
%   .thrs       - [nNodes x nTrees] threshold corresponding to each fid
%   .fids       - [nNodes x nTrees] feature ids for each node
%   .child      - [nNodes x nTrees] index of child for each node
%   .count      - [nNodes x nTrees] number of data points at each node
%   .depth      - [nNodes x nTrees] depth of each node
%   .eBins      - data structure for storing all node edge maps
%   .eBnds      - data structure for storing all node edge maps
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

% get default parameters
dfs={'imWidth',24, 'gtWidth',12, 'nPos',7e5, 'nNeg',7e5, 'nImgs',40, ...
  'nTrees',16, 'fracFtrs',1/8, 'minCount',1, 'minChild',8, ...
  'maxDepth',64, 'discretize','kmeans', 'nSamples',600, 'nClasses',2, ...
  'split','entropy', 'nOrients',0, 'grdSmooth',2, 'chnSmooth',2, ...
  'simSmooth',4, 'normRad',4, 'shrink',2, 'nCells',4, 'stride',2, ...
  'multiscale',0, 'sharpen',2, 'nTreesEval',10, 'nThreads',8, 'nms',0, ...
  'seed',1, 'useParfor',0, 'ctmaxval',5000, 'modelDir','models/', ...
  'nLandmarks',0, 'modelFnm','ctmodel', 'regSplit','mse', ...
  'imageDir','ct_training_data/', 'nodeSelectProb',0.0, 'stageId',0, ...
  'shpWidth',[], 'shpDepth',[], 'shpSmooth',1, 'shpShrink',4, 'nPosePar',0};

opts = getPrmDflt(varargin,dfs,1);
if(nargin==0), model=opts; return; end                                     

% if forest exists load it and return
cd(fileparts(mfilename('fullpath')));
forestDir = [opts.modelDir '/forest/'];
forestFn = [forestDir opts.modelFnm];
if(exist([forestFn '.mat'], 'file'))
  load([forestFn '.mat']); return; end

% compute constants and store in opts
nTrees=opts.nTrees; nCells=opts.nCells; shrink=opts.shrink;
opts.nPos=round(opts.nPos); opts.nNeg=round(opts.nNeg);
opts.nTreesEval=min(opts.nTreesEval,nTrees);
opts.stride=max(opts.stride,shrink);
imWidth=opts.imWidth; gtWidth=opts.gtWidth;
imWidth=round(max(gtWidth,imWidth)/shrink/2)*shrink*2;
opts.imWidth=imWidth; opts.gtWidth=gtWidth;
nChnsGrad=2+opts.nOrients; nChnsColor=1; 
nChns = nChnsGrad+nChnsColor; opts.nChns = nChns;
          
stageId=opts.stageId; opts.nShpFtrs=0; shpShrink=opts.shpShrink;
shpWidth=opts.shpWidth; shpDepth=opts.shpDepth;
shpWidth=round(shpWidth/shpShrink/2)*shpShrink*2;
shpDepth=round(shpDepth/shpShrink/2)*shpShrink*2;
opts.shpWidth=shpWidth; opts.shpDepth=shpDepth;
if(stageId),  opts.nShpFtrs = shpWidth*shpWidth*shpDepth/shpShrink/shpShrink/shpShrink; end
if(~stageId), opts.nPosePar = 0.0; end;

opts.nChnFtrs = imWidth*imWidth*imWidth*nChns/shrink/shrink/shrink;        
opts.nSimFtrs = (nCells*nCells*nCells)*(nCells*nCells*nCells-1)/2*nChns;   
opts.nTotFtrs = opts.nChnFtrs + opts.nSimFtrs + opts.nShpFtrs; disp(opts);

% generate stream for reproducibility of model
stream=RandStream('mrg32k3a','Seed',opts.seed);

% train nTrees random trees (can be trained with parfor if enough memory)
if(opts.useParfor), parfor i=1:nTrees, trainTree(opts,stream,i); end
else for i=1:nTrees, trainTree(opts,stream,i); end; end

% merge trees and save model
model = mergeTrees( opts );
if(~exist(forestDir,'dir')), mkdir(forestDir); end
save([forestFn '.mat'], 'model', '-v7.3');

end

function model = mergeTrees( opts )
% accumulate trees and merge into final model
nTrees=opts.nTrees; gtWidth=opts.gtWidth; 
numLm=opts.nLandmarks; nPosePar=opts.nPosePar;
treeFn = [opts.modelDir '/tree/' opts.modelFnm '_tree'];
for i=1:nTrees
  t=load([treeFn int2str2(i,3) '.mat'],'tree'); t=t.tree;
  if(i==1), trees=t(ones(1,nTrees)); else trees(i)=t; end
end
nNodes=0; for i=1:nTrees, nNodes=max(nNodes,size(trees(i).fids,1)); end

% merge all fields of all trees
model.opts=opts; Z=zeros(nNodes,nTrees,'uint32');
model.thrs     =zeros(nNodes,nTrees,'single');
model.gains    =zeros(nNodes,nTrees,'single');
model.meanOff  =zeros(numLm*3,nNodes,nTrees,'single');
model.covOff   =zeros(numLm*numLm*3*3,nNodes,nTrees,'single');
model.meanPose =zeros(nPosePar,nNodes,nTrees,'single');
model.covPose  =zeros(nPosePar*nPosePar,nNodes,nTrees,'single');
model.splitType=zeros(nNodes,nTrees,'uint8');
model.fids=Z; model.child=Z; model.count=Z; model.depth=Z;
model.segs=zeros(gtWidth,gtWidth,gtWidth,nNodes,nTrees,'uint8');
for i=1:nTrees, tree=trees(i); nNodes1=size(tree.fids,1);
  model.fids(1:nNodes1,i) = tree.fids;
  model.thrs(1:nNodes1,i) = tree.thrs;
  model.gains(1:nNodes1,i) = tree.gains;
  model.child(1:nNodes1,i) = tree.child;
  model.count(1:nNodes1,i) = tree.count;
  model.depth(1:nNodes1,i) = tree.depth;
  model.splitType(1:nNodes1,i) = tree.splitType;
  model.segs(:,:,:,1:nNodes1,i) = tree.hs-1;
  model.meanOff(:,1:nNodes1,i) = tree.meanOff;
  model.covOff(:,1:nNodes1,i) = tree.covOff;
  model.meanPose(:,1:nNodes1,i) = tree.meanPose;
  model.covPose(:,1:nNodes1,i) = tree.covPose;
end

% remove very small segments (<=5 pixels)
segs=model.segs; nSegs=squeeze(max(max(max(segs))))+1;
parfor i=1:nTrees*nNodes, m=nSegs(i);
  if(m==1), continue; end; S=segs(:,:,:,i); del=0;
  for j=1:m, Sj=(S==j-1); if(nnz(Sj)>5), continue; end
    S(Sj)=median(single(S(convTri3D(single(Sj),1)>0))); del=1; end
  if(del), [~,~,S]=unique(S); S=reshape(S-1,gtWidth,gtWidth,gtWidth);
    segs(:,:,:,i)=S; nSegs(i)=max(S(:))+1; end
end
model.segs=segs; model.nSegs=nSegs;

% store compact representations of sparse binary edge patches
nBnds=opts.sharpen+1; eBins=cell(nTrees*nNodes,nBnds);
eBnds=zeros(nNodes*nTrees,nBnds);
parfor i=1:nTrees*nNodes
  if(model.child(i) || model.nSegs(i)==1), continue; end %#ok<PFBNS>
  E=canny3D(single(model.segs(:,:,:,i)))>.01; E0=0;
  for j=1:nBnds, eBins{i,j}=uint16(find(E & ~E0)'-1); E0=E;
    eBnds(i,j)=length(eBins{i,j}); E=convTri3D(single(E),1)>.01; end
end
eBins=eBins'; model.eBins=[eBins{:}]';
eBnds=eBnds'; model.eBnds=uint32([0; cumsum(eBnds(:))]);
end

function trainTree( opts, stream, treeInd )
% Train a single tree in forest model.

% location of ground truth
trnImgDir = [opts.imageDir '/images/'];
trnGtDir  = [opts.imageDir '/groundtruth/'];
trnBouDir = [opts.imageDir '/boundary/'];
trnLMDir  = [opts.imageDir '/landmarks/'];
trnDofDir = [opts.imageDir '/dofs/'];
trnPemDir = [opts.imageDir '/pems/'];

imgIds=dir(trnImgDir); imgIds=imgIds([imgIds.bytes]>0); imgIds={imgIds.name};
extstr = numel(imgIds{1})-cell2array(regexp(imgIds{1},{'nii'})); 
ext = imgIds{1}(end-extstr:end);nImgs=length(imgIds); 
for i=1:nImgs, imgIds{i}=imgIds{i}(1:end-(extstr+2)); end

% check the number of landmarks 
fid = fopen([trnLMDir imgIds{1} '.txt']); lm=textscan(fid,'%f%f%f','delimiter','\n');
assert(opts.nLandmarks<=numel(lm{1})); numLm=opts.nLandmarks; avaLm=numel(lm{1}); 
fclose(fid); clear lm; 

% extract commonly used options
imWidth=opts.imWidth; imRadius=imWidth/2;
gtWidth=opts.gtWidth; gtRadius=gtWidth/2;
nChns=opts.nChns; nTotFtrs=opts.nTotFtrs;
nChnFtrs=opts.nChnFtrs; nSimFtrs=opts.nSimFtrs; nShpFtrs=opts.nShpFtrs;
nPos=opts.nPos; nNeg=opts.nNeg; shrink=opts.shrink;
nPosePar=opts.nPosePar;  stageId=opts.stageId;

% finalize setup
treeDir = [opts.modelDir '/tree/'];
treeFn = [treeDir opts.modelFnm '_tree'];
if(exist([treeFn int2str2(treeInd,3) '.mat'],'file'))
  fprintf('Reusing tree %d of %d\n',treeInd,opts.nTrees); return; end
fprintf('\n-------------------------------------------\n');
fprintf('Training tree %d of %d\n',treeInd,opts.nTrees); tStart=clock;

% set global stream to stream with given substream (will undo at end)
streamOrig = RandStream.getGlobalStream();
set(stream,'Substream',treeInd);
RandStream.setGlobalStream( stream );

% collect positive and negative patches and compute features
imgIds    = imgIds(randperm(nImgs,min(nImgs,opts.nImgs)));
k         = nPos+nNeg; nImgs=min(nImgs,opts.nImgs);
fids      = sort(randperm(nTotFtrs,round(nTotFtrs*opts.fracFtrs)));
ftrs      = zeros(k,length(fids),'single');
fWts      = zeros(1,length(fids),'single');
labels    = zeros(gtWidth,gtWidth,gtWidth,k,'uint8'); 
offsets   = zeros(3,numLm,k,'single'); 
posepars  = zeros(nPosePar,k,'single'); k = 0;
poseOrd   = {'rz','rx','ry'};
tid       = ticStatus('Collecting data',30,1);
spacing   = cell(nImgs,1);

for i = 1:nImgs
  
  % Load image, segmentation and boundaries
  I   = load_untouch_nii([trnImgDir imgIds{i} '.' ext]); I=I.img; siz=size(I);
  gt  = load_untouch_nii([trnGtDir  imgIds{i} '.' ext]); spacing{i}=gt.hdr.dime.pixdim(2:4); gt=gt.img;           
  bou = load_untouch_nii([trnBouDir imgIds{i} '.' ext]); bou=(bou.img)>0;  
  fid = fopen([trnLMDir imgIds{i} '.txt']); lm=textscan(fid,'%f%f%f','HeaderLines',avaLm-numLm,'delimiter','\n'); fclose(fid);
  I   = single(I)/single(opts.ctmaxval);                                      
  gt  = uint8(gt);
  if (ndims(I)~=3), error('trainTree:: currently supports only 3D images'); end;
  
  % Read the information for the hierarchical forest 
  if(stageId)
    pem = load_untouch_nii([trnPemDir imgIds{i} '.' ext]); pem=pem.img;
    vtk = vtk2Mat         ([trnPemDir imgIds{i} '_lm.vtk']);
    w2i = world2ImageMat  ([trnPemDir imgIds{i} '.' ext]);
    dof = readDofParMex   ([trnDofDir imgIds{i} '.dof.gz']);
  else
    pem=[];vtk=[];w2i=[];dof=[];
  end
  
  % Generate the boundary and segmentation structure
  gt = struct('Segmentation',gt,'Boundaries',bou,'Landmarks',{lm},'Dofs',dof); 
  clear bou; clear lm; clear dof; 
      
  % Perform padding and compute the channels
  p=mod(4-mod(siz(1:3),4),4);                                                % Changed - oo2113
  if(any(p)), I=padarray(I,p,'symmetric','post'); end                        % Changed - oo2113
  [chnsReg,chnsSim] = edgesChns_3D(I,opts);                                  % Changed - oo2113
  [chnsShp] = edgesShp_3D(pem,vtk,w2i,opts);
  
  % Sample positive and negative locations
  xyz=[]; k1=0; B=false(siz(1),siz(2),siz(3));                               % Changed - oo2113  
  B(shrink:shrink:end,shrink:shrink:end,shrink:shrink:end)=1;                % Changed - oo2113 
  B([1:imRadius end-imRadius:end],:,:)=0;
  B(:,[1:imRadius end-imRadius:end],:)=0;
  B(:,:,[1:imRadius end-imRadius:end])=0;

  % Pick the positive and negative sample indices
  M=gt.Boundaries; M(bwdist(M)<gtRadius)=1;
  [y,x,z]=ind2sub(siz,find(M.*B)); k2=min(length(y),ceil(nPos/nImgs));
  rp=randperm(length(y),k2); y=y(rp); x=x(rp); z=z(rp);
  xyz=[xyz; x y z ones(k2,1)*1]; k1=k1+k2;  p_t=k2; %#ok<AGROW> 
  [y,x,z]=ind2sub(siz,find(~M.*B)); k2=min(length(y),ceil(nNeg/nImgs));
  rp=randperm(length(y),k2); y=y(rp); x=x(rp); z=z(rp);
  xyz=[xyz; x y z zeros(k2,1)*1]; k1=k1+k2;  n_t=k2; %#ok<AGROW>
  if(k1>size(ftrs,1)-k), k1=size(ftrs,1)-k; xyz=xyz(1:k1,:); end
  display(sprintf('positive sample ratio is: %2.2f\n',p_t/(n_t+p_t)));
  
  % Crop patches and ground truth labels
  psReg=zeros(imWidth/shrink,imWidth/shrink,imWidth/shrink,nChns,k1,'single');
  psShp=cell(k1,1);
  lbls =zeros(gtWidth,gtWidth,gtWidth,k1,'uint8');
  offs =zeros(3,numLm,k1,'single');
  pose =zeros(nPosePar,k1,'single');
  [c_g,r_g,s_g]=meshgrid(0:siz(2)-1,0:siz(1)-1,0:siz(3)-1);  
  psSim=psReg; ri=imRadius/shrink; rg=gtRadius;
  
  for j=1:k1, xyz1=xyz(j,:); xyz2=xyz1/shrink;
    psReg(:,:,:,:,j)=chnsReg(xyz2(2)-ri+1:xyz2(2)+ri,xyz2(1)-ri+1:xyz2(1)+ri,xyz2(3)-ri+1:xyz2(3)+ri,:);
    psSim(:,:,:,:,j)=chnsSim(xyz2(2)-ri+1:xyz2(2)+ri,xyz2(1)-ri+1:xyz2(1)+ri,xyz2(3)-ri+1:xyz2(3)+ri,:);
    psShp{j} = chnsShp;
    t=gt.Segmentation(xyz1(2)-rg+1:xyz1(2)+rg,xyz1(1)-rg+1:xyz1(1)+rg,xyz1(3)-rg+1:xyz1(3)+rg);
    
    if (xyz1(4) == 1)
      offs(1,:,j)= gt.Landmarks{1}-r_g(xyz1(2),xyz1(1),xyz1(3));
      offs(2,:,j)= gt.Landmarks{2}-c_g(xyz1(2),xyz1(1),xyz1(3));
      offs(3,:,j)= gt.Landmarks{3}-s_g(xyz1(2),xyz1(1),xyz1(3));
    else
      offs(:,:,j)= zeros(3,numLm,'single');            
    end
    
    if(~isempty(gt.Dofs)), for p=1:nPosePar, pose(p,j) = getfield(gt.Dofs,poseOrd{p}); end; end%#ok<GFLD>
            
    if(all(t(:)==t(1))), lbls(:,:,:,j)=1; else [~,~,t]=unique(t);
      lbls(:,:,:,j)=reshape(t,gtWidth,gtWidth,gtWidth); end
  end  
  clear chnsReg; clear chnsSim; clear chnsShp; clear x_g; clear y_g; clear z_g;
  if(0), figure(1); montage2(squeeze(psReg(:,:,1,:))); drawnow; end
  if(0), figure(2); montage2(lbls(:,:,:)); drawnow; end
  
  % Compute features and store
  ftrs1=[reshape(psReg,[],k1)' stComputeSimFtrs(psSim,opts) reshape(cell2array(psShp),[],k1)'];
  fWts1=[ ones(nChnFtrs,1)*0.5; ones(nSimFtrs,1)*0.5; ones(nShpFtrs,1)*0.0 ];
  assert(size(ftrs1,2)==nTotFtrs); assert(isa(ftrs1,'single'));
  clear psSim; clear psReg; clear psShp;
  
  ftrs(k+1:k+k1,:)=ftrs1(:,fids);  clear ftrs1;
  fWts(:)=fWts1(fids);             clear fWts1;
  labels(:,:,:,k+1:k+k1)=lbls;     clear lbls;
  offsets(:,:,k+1:k+k1) =offs;     clear offs;
  posepars(:,k+1:k+k1)  =pose;     clear pose;
  k=k+k1; if(k==size(ftrs,1)), tocStatus(tid,1); break; end
  tocStatus(tid,i/nImgs);
end
if(k<size(ftrs,1)), ftrs=ftrs(1:k,:); labels=labels(:,:,:,1:k); offsets=offsets(:,:,1:k); posepars=posepars(:,1:k); end
spacing=cell2mat(spacing); assert(isrow(unique(spacing,'rows'))); rWts=repmat(spacing(1,:)',1,numLm); rWts=transpose(rWts(:));

% Train structured edge classifier (random decision tree)
pTree=struct('minCount',opts.minCount, 'minChild',opts.minChild, ...
  'maxDepth',opts.maxDepth, 'H',opts.nClasses, 'split',opts.split, ...
  'regSplit',opts.regSplit, 'nodeSelectProb',opts.nodeSelectProb, ...
  'rWts',rWts, 'fWts',fWts);

t=labels;   labels  =cell(k,1);                  for i=1:k, labels{i}=t(:,:,:,i); end; clear t;
t=offsets;  offsets =zeros(k,numLm*3,'single');  for i=1:k, t2=t(:,:,i); offsets(i,:)=t2(:); end; clear t; clear t2;
t=posepars; posepars=zeros(k,nPosePar,'single'); for i=1:k, posepars(i,:)=t(:,i); end; clear t;

% Pass the discretize function to the forest traning & train the forest
pTree.discretize=@(hs,H) discretize(hs,H,opts.nSamples,opts.discretize);
tree=hierHoughForestTrain(ftrs,labels,offsets,posepars,pTree);

% Save the tree & Correct selected feature ids 
tree.hs=cell2array(tree.hs);
tree.fids(tree.child>0) = fids(tree.fids(tree.child>0)+1)-1;
if(~exist(treeDir,'dir')), mkdir(treeDir); end
save([treeFn int2str2(treeInd,3) '.mat'],'tree'); e=etime(clock,tStart);
fprintf('Training of tree %d complete (time=%.1fs).\n',treeInd,e);
RandStream.setGlobalStream( streamOrig );

end

function ftrs = stComputeSimFtrs( chns, opts )
% Compute self-similarity features (order must be compatible w mex file).
w=opts.imWidth/opts.shrink; n=opts.nCells; if(n==0), ftrs=[]; return; end
nSimFtrs=opts.nSimFtrs; nChns=opts.nChns; m=size(chns,5);
inds=round(w/n/2); inds=round((1:n)*(w+2*inds-1)/(n+1)-inds+1); 
chns=reshape(chns(inds,inds,inds,:,:),n*n*n,nChns,m);
ftrs=zeros(nSimFtrs/nChns,nChns,m,'single');

k=0; for i=1:n*n*n-1, k1=n*n*n-i; i1=ones(1,k1)*i;
  ftrs(k+1:k+k1,:,:)=chns(i1,:,:)-chns((1:k1)+i,:,:); k=k+k1; end

ftrs = reshape(ftrs,nSimFtrs,m)';
end

function [hs,segs] = discretize( segs, nClasses, nSamples, type )

% Convert a set of segmentations into a set of labels in [1,nClasses].
persistent cache; w=size(segs{1},1); assert(size(segs{1},2)==w);  assert(size(segs{1},3)==w);
if(~isempty(cache) && cache{1}==w), [~,is1,is2]=deal(cache{:}); else
  % compute all possible lookup inds for w x w patches
  is=1:w^6; is1=floor((is-1)/w/w/w); is2=is-is1*w*w*w; is1=is1+1;
  kp=is2>is1; is1=is1(kp); is2=is2(kp); cache={w,is1,is2};
end

% compute n binary codes zs of length nSamples
nSamples=min(nSamples,length(is1)); kp=randperm(length(is1),nSamples);
n=length(segs); is1=is1(kp); is2=is2(kp); zs=false(n,nSamples);
for i=1:n, zs(i,:)=segs{i}(is1)==segs{i}(is2); end
zs=bsxfun(@minus,zs,sum(zs,1)/n); zs=zs(:,any(zs,1));
if(isempty(zs)), hs=ones(n,1,'uint32'); segs=segs{1}; return; end

% find most representative segs (closest to mean)
[~,ind]=min(sum(zs.*zs,2)); segs=segs{ind};

% apply PCA to reduce dimensionality of zs
U=pca(zs'); d=min(5,size(U,2)); zs=zs*U(:,1:d);

% discretize zs by clustering or discretizing pca dimensions
d=min(d,floor(log2(nClasses))); hs=zeros(n,1);
for i=1:d, hs=hs+(zs(:,i)<0)*2^(i-1); end
[~,~,hs]=unique(hs); hs=uint32(hs);
if(strcmpi(type,'kmeans'))
  nClasses1=max(hs); C=zs(1:nClasses1,:);
  for i=1:nClasses1, C(i,:)=mean(zs(hs==i,:),1); end
  hs=uint32(kmeans2(zs,nClasses,'C0',C,'nIter',1));
end

% optionally display different types of hs
for i=1:0, figure(i); montage2(cell2array(segs(hs==i))); end

end
