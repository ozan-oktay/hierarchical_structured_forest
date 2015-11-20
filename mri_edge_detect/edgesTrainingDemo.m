%% Demo for Structured Edge Detector (please see readme.txt first) %%%%
%%%% add the path
function edgesTrainingDemo()

restoredefaultpath;setenv('LD_LIBRARY_PATH','');
currentpath=pwd(); parsedpath=strsplit(currentpath,'/'); rootpath=strjoin(parsedpath(1:end-1),'/');
addpath(rootpath); addpath(genpath(strcat(rootpath,'/toolbox')));
%compileMex();  

%%%% set opts for training (see edgesTrain_3D.m) %%%%
   opts.imWidth= 28;                    % [28]  width of image patches    
   opts.gtWidth= 12;                    % [12]  width of ground truth patches     
      opts.nPos= 1.2e6;                 % [1e6] number of positive patches per tree
      opts.nNeg= 1.2e6;                 % [1e6] number of negative patches per tree
     opts.nImgs= Inf;                   % [75]  maximum number of images to use for training
    opts.nTrees= 8;                     % [16]  number of trees in forest to train
  opts.fracFtrs= 0.125;                 % [1/8] fraction of features to use to train each tree      
  opts.minCount= 1;                     % [1]   minimum number of data points to allow split
  opts.minChild= 8;                     % [8]   minimum number of data points allowed at child nodes
  opts.maxDepth= 64;                    % [64]  maximum depth of tree (default value 64)
opts.discretize= 'kmeans';              % ['kmeans'] options include 'pca' and 'kmeans'
  opts.nSamples= 600;                   % [600] number of samples for clustering structured labels 
  opts.nClasses= 2;                     % [2]   number of classes (clusters) for binary splits
     opts.split= 'entropy';             % ['entropy'] options include 'gini', 'entropy' and 'twoing'
  opts.nOrients= 0;                     % [0]   number of orientations per gradient scale   
 opts.grdSmooth= 2;                     % [2]   radius for image gradient smoothing (using convTri) 
 opts.chnSmooth= 2;                     % [2]   radius for reg channel smoothing (using convTri)
 opts.simSmooth= 4;                     % [4]   radius for sim channel smoothing (using convTri)
   opts.normRad= 4;                     % [4]   gradient normalization radius (see gradientMag)
    opts.shrink= 2;                     % [2]   amount to shrink channels
    opts.nCells= 4;                     % [4]   number of self similarity cells                       
    opts.stride= 2;                     % [2]   stride at which to compute edges
opts.multiscale= 1;                     % [0]   if true run multiscale edge detector
   opts.sharpen= 2;                     % [2]   sharpening amount (can only decrease after training)
opts.nTreesEval= 8;                     % [10]  number of trees to evaluate per location
  opts.nThreads= 8;                     % [8]   number of threads for evaluation of trees
       opts.nms= 0;                     % [0]   if true apply non-maximum suppression to edges
      opts.seed= 1;                     % [1]   seed for random stream (for reproducibility)
 opts.useParfor= 1;                     % [0]   if true train trees in parallel (memory intensive)
  opts.modelDir= 'models/';             % ['models/'] target directory for storing models
  opts.modelFnm= 'new_prn016_shape_large';    % ['ctmodel'] model filename
  opts.imageDir= 'mritrainingdata_sec_large/';% ['mritrainingdata_sec/'] location of image dataset     
  opts.ctmaxval= 1024;                  % [1024] maximum allowed intensity value - for linear scaling. 
  
opts.nLandmarks= 6;                     % [false] if true train trees with both classification and regression nodes
  opts.regSplit= 'mse';                 % ['mse'] regression criterion
  opts.nodeProb= [0.160,0.420,0.420];   % probability of selecting a regression split in training
   opts.stageId= 1;                     % Framework stage identifier - the zeroth level extracts pems/landmarks - the first level outputs also pose/scale in addition 
  opts.shpWidth= 102;                   % Width of the window for PEM Shape features
  opts.shpDepth= 18;                    % Depth of the window for PEM Shape features     
 opts.shpSmooth= 3;                     % Smoothing sigma par for PEM Shape features
 opts.shpShrink= 3;                     % Shrink factor for PEM Shape features
 
   opts.nPosePar= 2;                    % Number of pose parameters used in tree training (grountruth information)
opts.nShpOrients= 8;                    % Number of orientations in HoG histograms
opts.nShpBinSize= 12;                   % Number of pixels in each dim for histogram binning
  
%%%% train edge detector (~150m/20Gb per tree, proportional to nPos/nNeg) %%%%
if (~isempty(gcp('nocreate'))), poolobj = gcp('nocreate'); delete(poolobj); end; 
if (opts.useParfor), pool=parpool(3); end;
tic, edgesTrain_3D(opts); toc; % will load model if already trained
if (~isempty(gcp('nocreate'))), poolobj = gcp('nocreate'); delete(poolobj); end; 
