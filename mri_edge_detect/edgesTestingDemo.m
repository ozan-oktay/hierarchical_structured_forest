%%%% Demo for Structured Edge Detector (please see readme.txt first).

function edgesTestingDemo (modelname,imagename,savename)
  %% set opts for testing and training (see edgesTrain_3D.m) %%%%
  restoredefaultpath;setenv('LD_LIBRARY_PATH',''); 
  workdirectory = '/vol/biomedic/users/oo2113/str_hier_forest_mri'; addpath(workdirectory);
  addpath(genpath(horzcat(workdirectory,'/toolbox')));
  compileMex();  
  
  opts=edgesTrain_3D();                    % default options (good settings)
  opts.modelDir   = 'models/';             % model will be in models/forest
  opts.modelFnm   = 'mriSecond_hier_X';    % model name
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
  opts.imagename          = horzcat(workdirectory,'/mritestingdata/images/14DB02502_ed0_3D.nii.gz');
  opts.patname            = strsplit(opts.imagename,'/'); opts.patname=opts.patname{end};
  opts.savename           = horzcat(workdirectory,sprintf('/tmp/n_pem_%s.nii.gz',opts.patname));
  opts.houghimagesavename = horzcat(workdirectory,sprintf('/tmp/n_hou_%s.nii.gz',opts.patname));
  opts.vtksavename        = horzcat(workdirectory,sprintf('/tmp/n_lm_%s.vtk',opts.patname));
  opts.posetxtsavename    = horzcat(workdirectory,sprintf('/tmp/n_pose_%s.txt',opts.patname));
  opts.tmpimagename       = strcat(tempname('./tmp'),'.nii.gz');
  opts.referenceimagename = horzcat(workdirectory,'/',model.opts.imageDir,'/reference/ref_image.nii.gz');
  if (nargin >= 2), opts.imagename = imagename; end;
  if (nargin >= 3), opts.savename = savename;  filebasename=strsplit(savename,'.'); filebasename=filebasename{1};
    opts.houghimagesavename=strcat(filebasename,'_lm.nii.gz'); opts.vtksavename=strcat(filebasename,'_lm.vtk');opts.posetxtsavename=strcat(filebasename,'_pose.txt'); end;
  
  %% if the classifier is in the second stage find the pem image %%%%
  if (model.opts.stageId == 1)
    tmp1              = strsplit(opts.imagename,'/');patName = tmp1{end};tmp2=strsplit(patName,'.');
    opts.pemimagename = strcat(strjoin(tmp1(1:end-2),'/'),'/pems/',patName);
    opts.vtkname      = strcat(strjoin(tmp1(1:end-2),'/'),'/pems/',tmp2{1},'_lm.vtk');
    otherFtrs.pem     = load_untouch_nii(opts.pemimagename); assert(all(otherFtrs.pem.hdr.dime.pixdim(2:4)==[1.25,1.25,2.00])); otherFtrs.pem = otherFtrs.pem.img; 
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
  save_untouch_nii(Imgstr,opts.savename);
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
  fprintf (sprintf('Estimated pose values (cov): %3.3f %3.3f\n' ,pVal(2),pVal(4)));
  fileID = fopen(opts.posetxtsavename,'w');fprintf(fileID,'mean val: %3.3f %3.3f\n',pVal(1),pVal(3));fprintf(fileID,'cov val: %3.3f %3.3f\n',pVal(2),pVal(4));fclose(fileID);
    
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
