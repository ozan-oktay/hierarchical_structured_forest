function analyze_model()
close all;
model=[];
colors=['-*b';'-*g';'-*r';'-*k';'-*c'];
fprintf('Loading the model... \n');
currentpath=pwd(); parsedpath=strsplit(currentpath,'/'); rootpath=strjoin(parsedpath(1:end-1),'/');
addpath(rootpath); modelName=strcat(rootpath,'/models/forest/mriSecond_hier_AFFT.mat');
load(modelName); addpath(genpath(strcat(rootpath,'/toolbox')));
fprintf('Model is loaded.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY THE PROXIMITY MAP OF TRAINING IMAGES
if (~isempty(model.dataInf{1}.affMat))

  nTrees   = model.opts.nTrees;
  for ii=2:nTrees, model.dataInf{1}.affMat=model.dataInf{1}.affMat+model.dataInf{ii}.affMat;end
  affMat   = model.dataInf{1}.affMat / nTrees;
  nImgs    = size(affMat,1); 
  rotmat   = strcat(rootpath,'/mritrainingdata_sec/dofs/rotation.mat'); load(rotmat,'filename','z_rot');
  filename = cellstr(filename); %#ok<NODEF>
  rotval   = cell(nImgs,1);
  
  %%%%%%%% FIND THE ROTATION VALUES OF THE TRAINING DATASET %%%%%%%%
  for mm=1:nImgs
       t_name=model.dataInf{1}.imgIds{mm}; res=strfind(filename,t_name); nn=find(~cellfun(@isempty,res)); rotval{mm}=z_rot(nn);     
  end

  %%%%%%%% CONVERT IT TO A DISSIMILARITY MATRIX %%%%%%%%
  rotval = squeeze(cell2array(rotval));
  affMat = 1./(affMat+1e-15); affMat(eye(nImgs)>0)=0.0;             
  [Y,E] = cmdscale(affMat);

  figure(1);
  scatter(Y(:,1),Y(:,2),30,rotval,'filled');colormap('jet'); colorbar; grid on;
  xlabel('Dimension 1'); ylabel('Dimension 2'); title('Proximity Plot of Images VS Rotation Info');
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY THE DISTRIBUTION OF SPLIT TYPES VS TREE DEPTH
nDepth        = min(max(model.depth));
nSplitTypes   = max(max(model.splitType));
splitTypeDist = zeros(nDepth,nSplitTypes);
for dId=1:nDepth
  indices    = (model.depth==dId);
  splitTypes = model.splitType(indices);
  nElem      = histc(splitTypes,1:nSplitTypes);
  for sId=1:nSplitTypes
    splitTypeDist(dId,sId)=nElem(sId)/sum(nElem);
  end
end

figure(2); 
for sId=1:nSplitTypes, plot(1:nDepth,splitTypeDist(:,sId),colors(sId,:),'LineWidth',2); hold on; end
grid on; h_legend=legend('Classification Node','Offset Regression Node','Rotation Regression Node','Location','NorthWest'); set(h_legend,'FontSize',14);
xlabel('Tree Depth'); ylabel('Perc of Nodes');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY THE USAGE OF FEATURE TYPES VS TREE DEPTH
nChnFtrs      = model.opts.nChnFtrs;
nSimFtrs      = model.opts.nSimFtrs;
nShpFtrs      = model.opts.nShpFtrs;
featureLim    = [0,nChnFtrs,nChnFtrs+nSimFtrs,nChnFtrs+nSimFtrs+nShpFtrs]+1;
nDepth        = min(max(model.depth));
nFeatureTypes = 3;
fidDist       = zeros(nDepth,nFeatureTypes);

for dId=1:nDepth
  indices  = (model.depth==dId);
  nodeFids = model.fids(indices);
    nElem  = histc(nodeFids,featureLim);
  for sId=1:nFeatureTypes
    fidDist(dId,sId)=nElem(sId)/sum(nElem);
  end
end

figure(3); 
for sId=1:nFeatureTypes, plot(1:nDepth,fidDist(:,sId),colors(sId,:),'LineWidth',2); hold on; end
grid on; h_legend=legend('Channel Features','SelfSimilarity Features','Shape Features','Location','NorthWest'); set(h_legend,'FontSize',14);
xlabel('Tree Depth'); ylabel('Perc of Nodes');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY IF THE TREE IS FULL OR NOT 
nDepth   = min(max(model.depth));
nodeDist = zeros(nDepth,2);
for dId=1:nDepth 
  expNumNodes     = pow2(double(dId));
  numTrees        = double(model.opts.nTrees);
  compNumNodes    = sum(model.depth(:)==dId)/numTrees;
  nodeDist(dId,1) = log2(compNumNodes);
  nodeDist(dId,2) = log2(expNumNodes);
end
figure(4); plot(1:nDepth,nodeDist(:,1),colors(1,:),'LineWidth',2); hold on;
           plot(1:nDepth,nodeDist(:,2),colors(2,:),'LineWidth',2); grid on; xlabel('Tree Depth'); ylabel('Log 2 of Number of Nodes at each level');
           h_legend=legend('Computed Num Nodes','Expected Num Nodes','Location','NorthWest'); set(h_legend,'FontSize',14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY HOW MUCH GAIN IS OBTAINED FROM EACH FEATURE
nDepth   = min(max(model.depth));
nChnFtrs      = model.opts.nChnFtrs;
nSimFtrs      = model.opts.nSimFtrs;
nShpFtrs      = model.opts.nShpFtrs;
nFeatureTypes = 3;
gainDist      = zeros(nDepth,nFeatureTypes);
featureLim    = [0,nChnFtrs,nChnFtrs+nSimFtrs,nChnFtrs+nSimFtrs+nShpFtrs]+1;

for dId=1:nDepth
  indices  = (model.depth==dId);
  fidsInd  = model.fids(indices);
  gainInd  = model.gains(indices);
  for sId=1:nFeatureTypes
    nodeIndices = (featureLim(sId)<=fidsInd) & (fidsInd<featureLim(sId+1));
    gainDist(dId,sId) = sum(gainInd(nodeIndices)) / numel(gainInd(nodeIndices));  
  end
end

figure(5); 
for sId=1:nFeatureTypes, plot(1:nDepth,gainDist(:,sId),colors(sId,:),'LineWidth',2); hold on; end
grid on; h_legend=legend('Channel Features','SelfSimilarity Features','Shape Features','Location','NorthWest'); set(h_legend,'FontSize',14);
xlabel('Tree Depth'); ylabel('Average Information Gain');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end