function analyze_model()
close all;
model=[];
colors=['-*b';'-*g';'-*r'];
inv_colors=['-*r';'-*g';'-*b'];
fprintf('Loading the model... \n');
currentpath=pwd(); parsedpath=strsplit(currentpath,'/'); rootpath=strjoin(parsedpath(1:end-1),'/');
addpath(rootpath); modelName=strcat(rootpath,'/models/forest/myforest_olddata.mat');
load(modelName); addpath(genpath(strcat(rootpath,'/toolbox')));
fprintf('Model is loaded.\n');
display(model.opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY THE PROXIMITY MAP OF TRAINING IMAGES
if (~isempty(model.dataInf{1}.affMat))

  % Build the affinity matrix  
  nTrees = model.opts.nTrees;
  for ii=2:nTrees, model.dataInf{1}.affMat=model.dataInf{1}.affMat+model.dataInf{ii}.affMat;end

  affMat   = model.dataInf{1}.affMat / nTrees;
  nImgs    = size(affMat,1); 
  rotmat   = strcat(rootpath,'/',model.opts.imageDir,'dofs/patientParam.mat'); load(rotmat,'filename','z_rot','scale'); 
  filename = cellstr(filename); %#ok<NODEF>
  rotval   = cell(nImgs,1);
  scaleval = cell(nImgs,1);
  
  %%%%%%%% FIND THE ROTATION VALUES OF THE TRAINING DATASET %%%%%%%%
  for mm=1:nImgs
       t_name=model.dataInf{1}.imgIds{mm}; res=strfind(filename,t_name); nn=find(~cellfun(@isempty,res)); rotval{mm}=z_rot(nn); scaleval{mm}=scale(nn);
  end

  %%%%%%%% CONVERT IT TO A DISSIMILARITY MATRIX %%%%%%%%
  rotval   = squeeze(cell2array(rotval));
  scaleval = squeeze(cell2array(scaleval)); 
  affMat   = 1./(affMat+1e-15); affMat(eye(nImgs)>0)=0.0;             
  [Y,E]    = leigs(affMat, 25, 3); display(sprintf('eigenvalues are: %f ',diag(E)));
  
  %%%% Plot The affinities %%%%
  figure(1);
  scatter(Y(:,1),Y(:,2),50,rotval,'filled');colormap('jet'); colorbar; grid on;
  xlabel('Dimension 1'); ylabel('Dimension 2'); title('Proximity Plot of Images VS Rotation Info');
  
  figure(7);
  scatter(Y(:,1),Y(:,3),50,scaleval,'filled');colormap('jet'); colorbar; grid on;
  xlabel('Dimension 1'); ylabel('Dimension 3'); title('Proximity Plot of Images VS Scale Info');
  
  figure(8);
  scatter(Y(:,2),Y(:,3),50,scaleval,'filled');colormap('jet'); colorbar; grid on;
  xlabel('Dimension 2'); ylabel('Dimension 3'); title('Proximity Plot of Images VS Scale Info');
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
grid on; h_legend=legend('Stratification Node','Classification Node','Regression Node','Location','NorthWest'); set(h_legend,'FontSize',14);
xlabel('Tree Depth'); ylabel('Perc of Nodes');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY THE USAGE OF FEATURE TYPES VS TREE DEPTH
%{
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
for sId=1:nFeatureTypes, plot(1:nDepth,fidDist(:,sId),inv_colors(sId,:),'LineWidth',2); hold on; end
grid on; h_legend=legend('Channel Features','SelfSimilarity Features','Shape Features','Location','NorthWest'); set(h_legend,'FontSize',14);
xlabel('Tree Depth'); ylabel('Perc of Nodes');
%}
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
% DISPLAY HOW MUCH AVERAGE GAIN IS OBTAINED FOR EACH SPLIT TYPE
%{
nDepth        = min(max(model.depth));
nSplitTypes   = max(max(model.splitType));
gainDist      = zeros(nDepth,nSplitTypes);

for dId=1:nDepth
  indices    = (model.depth==dId);
  splitTypes = model.splitType(indices);
  gainInd    = model.gains(indices);
  for sId=1:nSplitTypes
    nodeIndices       = (splitTypes==sId); 
    gainDist(dId,sId) = sum(gainInd(nodeIndices)) / numel(gainInd(nodeIndices));  
  end
end

figure(5); 
for sId=1:nSplitTypes, plot(1:nDepth,gainDist(:,sId),colors(sId,:),'LineWidth',2); hold on; end
grid on; h_legend=legend('Stratification Node','Classification Node','Regression Node','Location','NorthWest'); set(h_legend,'FontSize',14);
xlabel('Tree Depth'); ylabel('Average Information Gain');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOW MUCH AVERAGE GAIN IS OBTAINED FOR EACH INDIVIDUAL STRATIFICATION FOREST 
nDepth        = min(max(model.depth));
nSplitTypes   = model.opts.nPosePar;
gainDist      = zeros(nDepth,nSplitTypes+1);

for dId=1:nDepth
  indices    = (model.depth==dId);
  splitTypes = model.splitType(indices);
  for sId=1:nSplitTypes
    gainIndices       = model.gainsInd(sId,indices);
    nodeIndices       = (splitTypes==1); 
    gainDist(dId,sId) = sum(gainIndices(nodeIndices)) / numel(gainIndices(nodeIndices));  
  end
  gainIndices      = model.gains(indices);
  nodeIndices      = (splitTypes==1);
  gainDist(dId,end) = sum(gainIndices(nodeIndices)) / numel(gainIndices(nodeIndices));  
end

figure(6); 
for sId=1:nSplitTypes+1, plot(1:nDepth,1-gainDist(:,sId),colors(sId,:),'LineWidth',2); hold on; end
grid on; h_legend=legend('Rotation Part','Scale Part','Joint Stratification','Location','NorthWest'); set(h_legend,'FontSize',14);
xlabel('Tree Depth'); ylabel('Average Entropy');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end




% --- leigs function for Laplacian eigenmap.
% Written by Belkin & Niyogi, 2002.
function [E,V] = leigs(inpAff, PARAM, NE) 

  [Z,I] = sort ( inpAff,2);
  n     = size(inpAff,1);
  A     = zeros(n,n); 
  for i=1:n
    for j=2:PARAM+1
        A(i,I(i,j))= Z(i,j); 
        A(I(i,j),i)= Z(i,j); 
    end;    
  end;

  display('Adjacency Matrix is ready'); clear Z; clear I; clear dt;
  W = A; 

  [A_i, A_j, A_v] = find(A);  % disassemble the sparse matrix
  clear A; 
  t = mean(A_v)*2;

  for i = 1: size(A_i)  
      W(A_i(i), A_j(i)) = exp(1).^((-1/t)*(A_v(i)));
      %W(A_i(i), A_j(i))  = 1.0;
  end;
  D = sum(W(:,:),2);   

  L = spdiags(D,0,speye(size(W,1)))-W;
  display('Laplacian Graph is Ready')
  clear D; clear W; clear A_i; clear A_j; clear A_v;
  opts.tol = 1e-9;
  opts.issym=1; 
  opts.disp = 0; 
  [E,V] = eigs(L,NE+1,'sm',opts);
  V = diag(V); V = diag(V(NE:-1:1));
  E = E(:,NE:-1:1);
  
end