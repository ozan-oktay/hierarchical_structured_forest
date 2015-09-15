%%%% Demo for Structured Edge Detector (please see readme.txt first).

function edgesTestingDemo (modelname,imagename,savedir)
  %% set opts for testing and training (see edgesTrain_3D.m) %%%%
  restoredefaultpath;setenv('LD_LIBRARY_PATH',''); 
  currentpath=pwd(); parsedpath=strsplit(currentpath,'/'); workdirectory=strjoin(parsedpath(1:end-1),'/');
  addpath(workdirectory); addpath(genpath(strcat(workdirectory,'/toolbox')));
  %compileMex();  
  
  opts=edgesTrain_3D();                    % default options (good settings)
  opts.modelDir   = 'models/';             % model will be in models/forest
  opts.modelFnm   = 'mriFirst_hier_X2';    % model name
  if (nargin >= 1), opts.modelFnm  = modelname; end;
  
  model=edgesTrain_3D(opts);           	     % will load model if already trained
  model.opts.multiscale = 1;             	   % for top accuracy set multiscale=1
  model.opts.stride     = 1;                 % stride at which to compute edges
  model.opts.sharpen    = 1;       	         % for top speed set sharpen=0
  model.opts.nTreesEval = model.opts.nTrees; % for top speed set nTreesEval=1
  model.opts.nms        = -1;          	     % set to true to enable nms
  
  otherFtrs               = struct('pem',[], 'vtk',[]);     % pem image name - required for hier forest
  opts.pixel_type         = 'uint16';         % output image pixel type
  opts.pixel_spacing      = [1.25,1.25,2.00]; % intensity and pem pixel spacing in world coordinates
  opts.referenceimagename = horzcat(workdirectory,'/',model.opts.imageDir,'reference/ref_image.nii.gz');
  opts.imagename          = horzcat(workdirectory,'/mritestingdata/images/14DA01414_ed0_3D.nii.gz');
  opts.savedir            = horzcat(workdirectory,'/tmp');  

  if (nargin >= 2), opts.imagename = imagename; end;
  if (nargin >= 3), opts.savedir   = savedir;   end; mkdir(opts.savedir)
  [~,opts.patname,~]      =fileparts(opts.imagename); t=strsplit(opts.patname,'.nii'); opts.patname=t{1};
  opts.pemsavename        =strcat(opts.savedir,'/',opts.patname,'_pem.nii.gz'); 
  opts.houghimagesavename =strcat(opts.savedir,'/',opts.patname,'_lm.nii.gz'); 
  opts.vtksavename        =strcat(opts.savedir,'/',opts.patname,'_lm.vtk');
  opts.posetxtsavename    =strcat(opts.savedir,'/',opts.patname,'_pose.txt'); 
  opts.tmpimagename       =strcat(tempname('./tmp'),'.nii.gz');
  
  %% if the classifier is in the second stage find the pem image %%%%
  if (model.opts.stageId == 1)
    tmp1              = strsplit(opts.savedir,'/');    
    opts.pemimagename = strcat(strjoin(tmp1(1:end-1),'/'),'/pems/',opts.patname,'_pem.nii.gz');
    opts.vtkname      = strcat(strjoin(tmp1(1:end-1),'/'),'/pems/',opts.patname,'_lm.vtk');
    otherFtrs.pem     = load_untouch_nii(opts.pemimagename); assert(all((otherFtrs.pem.hdr.dime.pixdim(2:4)-[1.25,1.25,2.00])<1e-4)); otherFtrs.pem = otherFtrs.pem.img; 
    otherFtrs.w2i     = world2ImageMat  (opts.pemimagename);
    otherFtrs.vtk     = vtk2Mat         (opts.vtkname);
  end 
  
  %% detect the edges in the image %%%%
  funPreprocessing (opts);
  Imgstr                  = load_untouch_nii(opts.tmpimagename);
  Imgstr.img              = single(Imgstr.img);
  [outputImage,~,HV,pVal] = edgesDetect_3D(Imgstr.img/single(model.opts.ctmaxval),model,otherFtrs);
  
  %%%% save the results (.nii.gz) %%%%
  Imgstr.hdr.dime.datatype = int16(4);
  Imgstr.img = cast(1024 * outputImage, opts.pixel_type);
  save_untouch_nii(Imgstr,opts.pemsavename);
  delete(opts.tmpimagename);

  %%%% save the hough votes %%%%
  Imgstr.hdr.dime.datatype = int16(16);
  Imgstr.hdr.dime.dim(1)   = 4;
  Imgstr.hdr.dime.dim(5)   = model.opts.nLandmarks;
  Imgstr.img               = cast(HV, 'single');
  save_untouch_nii(Imgstr,opts.houghimagesavename);
  
  %%%% find the landmark locations from the hough vote maps 
  vtkParams=struct( 'method','max', 'savename',opts.vtksavename, 'refimagename',opts.imagename, 'spacing',Imgstr.hdr.dime.pixdim(2:4) );
  tstart=tic; houghVotes2Vtk (HV,vtkParams); fprintf (sprintf('time lm extraction: %f sec\n',toc(tstart)));
  
  % save the pose values 
  fprintf (sprintf('Estimated pose values (mean): %3.3f %3.3f\n',pVal(1),pVal(3)));
  fprintf (sprintf('Estimated pose values (std): %3.3f %3.3f\n' ,sqrt(pVal(2)),sqrt(pVal(4))));
  fileID = fopen(opts.posetxtsavename,'w');fprintf(fileID,'mean val: %3.3f %3.3f\n',pVal(1),pVal(3));fprintf(fileID,'std val: %3.3f %3.3f\n',sqrt(pVal(2)),sqrt(pVal(4)));fclose(fileID);
    
end

%% Preprocessing function for the input images / The same procedure is
%%%% required for the training data
function funPreprocessing (op)

%%%% Resample to fixed isotropic resolution 
binaryresample  ='irtk_binaries/linux/resample';
system(sprintf('%s %s %s -size %f ',binaryresample,op.imagename,op.tmpimagename, op.pixel_spacing)); 

%%%% Histogram Normalization ( not for CT images - other modalities )
binarynormalize ='irtk_binaries/linux/normalize';
system(sprintf('%s %s %s %s -piecewise',binarynormalize,op.referenceimagename,op.tmpimagename,op.tmpimagename)); 

%%%% Convert the images to NIFTI format (can be used to convert to NIFTI format)
binaryconvert   ='irtk_binaries/linux/convert';
system(sprintf('%s %s %s -ushort',binaryconvert,op.tmpimagename,op.tmpimagename)); 

end
