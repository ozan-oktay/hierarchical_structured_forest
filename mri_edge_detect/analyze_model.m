function analyze_model()
%close all;
colors=['-*b';'-*g';'-*r';'-*k';'-*c'];
fprintf('Loading the model... \n');
modelName='/vol/biomedic/users/oo2113/str_hier_forest_mri/models/forest/mriSecond_hier_E4.mat';
load(modelName);
fprintf('Model is loaded.\n');

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

figure(1); 
for sId=1:nSplitTypes, plot(1:nDepth,splitTypeDist(:,sId),colors(sId,:),'LineWidth',2); hold on; end
grid on; h_legend=legend('Classification Node','Offset Regression Node','Rotation Regression Node'); set(h_legend,'FontSize',14);
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

figure(2); 
for sId=1:nFeatureTypes, plot(1:nDepth,fidDist(:,sId),colors(sId,:),'LineWidth',2); hold on; end
grid on; h_legend=legend('Channel Features','SelfSimilarity Features','Shape Features'); set(h_legend,'FontSize',14);
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
figure(3); plot(1:nDepth,nodeDist(:,1),colors(1,:),'LineWidth',2); hold on;
           plot(1:nDepth,nodeDist(:,2),colors(2,:),'LineWidth',2); grid on; xlabel('Tree Depth'); ylabel('Log 2 of Number of Nodes at each level');
           h_legend=legend('Computed Num Nodes','Expected Num Nodes'); set(h_legend,'FontSize',14);
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
    gainDist(dId,sId) = sum(gainInd(nodeIndices));  
  end
end

figure(4); 
for sId=1:nFeatureTypes, plot(1:nDepth,gainDist(:,sId),colors(sId,:),'LineWidth',2); hold on; end
grid on; h_legend=legend('Channel Features','SelfSimilarity Features','Shape Features'); set(h_legend,'FontSize',14);
xlabel('Tree Depth'); ylabel('Information Gain');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

