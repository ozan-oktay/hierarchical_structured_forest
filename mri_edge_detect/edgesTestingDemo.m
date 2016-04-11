%%%% Demo for Structured Edge Detector (please see readme.txt first).

function edgesTestingDemo (modelname,imagename,savedir)
  %% set opts for testing and training (see edgesTrain_3D.m) %%%%
  restoredefaultpath;setenv('LD_LIBRARY_PATH',''); 
  currentpath=pwd(); parsedpath=strsplit(currentpath,'/'); workdirectory=strjoin(parsedpath(1:end-1),'/');
  addpath(workdirectory); addpath(genpath(strcat(workdirectory,'/toolbox')));
  %compileMex();  
  
  opts=edgesTrain_3D();                     % default options (good settings)
  opts.modelDir   = 'models/';              % model will be in models/forest
  opts.modelFnm   = 'myforest_firstStage';  % model name
  if (nargin >= 1), opts.modelFnm  = modelname; end;
  
  model=edgesTrain_3D(opts);           	     % will load model if already trained
  model.opts.multiscale = 0;             	   % for top accuracy set multiscale=1
  model.opts.stride     = 1;                 % stride at which to compute edges
  model.opts.sharpen    = 1;       	         % for top speed set sharpen=0
  model.opts.nTreesEval = model.opts.nTrees; % for top speed set nTreesEval=1
  model.opts.nms        = -1;          	     % set to true to enable nms
  
  otherFtrs               = struct('pem',[], 'vtk',[]);     % pem image name - required for hier forest
  opts.pixel_type         = 'uint16';         % output image pixel type
  opts.pixel_spacing      = [1.25,1.25,2.00]; % intensity and pem pixel spacing in world coordinates
  opts.referenceimagename = horzcat(workdirectory,'/',model.opts.imageDir,'/reference/ref_image.nii.gz');
  opts.referencecropname  = horzcat(workdirectory,'/',model.opts.imageDir,'/reference/ref_image_cropped.nii.gz');
  opts.imagename          = horzcat(workdirectory,'/mritestingdata/images/14DA01414_ed0_3D.nii.gz');
  opts.savedir            = horzcat(workdirectory,'/tmp/pems');  
  opts.stageId            = model.opts.stageId;

  if (nargin >= 2), opts.imagename = imagename; end;
  if (nargin >= 3), opts.savedir   = savedir;   end; mkdir(opts.savedir)
  [~,opts.patname,~]      =fileparts(opts.imagename); t=strsplit(opts.patname,'.nii'); opts.patname=t{1};
  opts.pemsavename        =strcat(opts.savedir,'/',opts.patname,'_pem.nii.gz'); 
  opts.houghimagesavename =strcat(opts.savedir,'/',opts.patname,'_lm.nii.gz'); 
  opts.vtksavename        =strcat(opts.savedir,'/',opts.patname,'_lm.vtk');
  opts.posetxtsavename    =strcat(opts.savedir,'/',opts.patname,'_pose.txt'); 
  opts.entropysavename    =strcat(opts.savedir,'/',opts.patname,'_entropy.txt'); 
  opts.tmpimagename       =strcat(tempname('./tmp'),'.nii.gz');
  
  %% if the classifier is in the second stage find the pem image %%%%
  if (opts.stageId == 1)
    tmp1              = strsplit(opts.savedir,'/');    
    opts.pemimagename = strcat(strjoin(tmp1(1:end-1),'/'),'/pems/',opts.patname,'_pem.nii.gz');
    opts.vtkname      = strcat(strjoin(tmp1(1:end-1),'/'),'/pems/',opts.patname,'_lm.vtk');
    otherFtrs.pem     = load_untouch_nii(opts.pemimagename); assert(all(abs(otherFtrs.pem.hdr.dime.pixdim(2:4)-opts.pixel_spacing)<1e-4)); otherFtrs.pem = otherFtrs.pem.img; 
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

  %%%% save the hough votes %%%%
  Imgstr.hdr.dime.datatype = int16(16);
  Imgstr.hdr.dime.dim(1)   = 4;
  Imgstr.hdr.dime.dim(5)   = model.opts.nLandmarks;
  Imgstr.img               = cast(HV, 'single');
  save_untouch_nii(Imgstr,opts.houghimagesavename);
  
  %%%% find the landmark locations from the hough vote maps 
  vtkParams=struct( 'method','max', 'savename',opts.vtksavename, 'refimagename',opts.tmpimagename, 'spacing',Imgstr.hdr.dime.pixdim(2:4) );
  tstart=tic; houghVotes2Vtk (HV,vtkParams); fprintf (sprintf('time lm extraction: %f sec\n',toc(tstart)));
  delete(opts.tmpimagename);

  % save the pose values 
  if (opts.stageId == 1)
      fileID = fopen(opts.posetxtsavename,'w');
      fprintf(fileID,'mean val: %3.3f %3.3f\n',pVal(1),pVal(3));
      fprintf(fileID,'std val: %3.3f %3.3f\n' ,sqrt(pVal(2)),sqrt(pVal(4)));
      fclose(fileID);
  end
  
  %% ESTIMATE THE CONFIDENCE VALUE %%
  %%% Google Search: measure agreement of hough votes /or/ consistency
%   conf_values=zeros(model.opts.nLandmarks,1);
%   max_prob_val=zeros(model.opts.nLandmarks,1);
%   for channelId=1:size(HV,4)
%     HV_chn  = HV(:,:,:,channelId);crop_size=30;
%     index = find(HV_chn == max(HV_chn(:))); [r_m,c_m,s_m] = ind2sub(size(HV_chn),index);
%     %r_b = (max(r_m-crop_size,1):min(r_m+crop_size,size(HV_chn,1)));
%     %c_b = (max(c_m-crop_size,1):min(c_m+crop_size,size(HV_chn,2)));
%     %s_b = (max(s_m-crop_size,1):min(s_m+crop_size,size(HV_chn,3))); HV_chn = HV_chn(r_b,c_b,s_b);
%     HV_chn = HV_chn / sum(HV_chn(:));    
%     siz = size(HV_chn);
%     [row,col,slc] = ndgrid(1:siz(1),1:siz(2),1:siz(3)); 
%     par      = [row(:),col(:),slc(:)]; numPar=size(par,2);
%     cov_mat  = zeros(numPar,numPar);
%     mean_mat = zeros(numPar,1);
%     % Mean Vector and Covariance Matrix Computation
%     for parId1=1:numPar, mean_mat(parId1) = sum( par(:,parId1) .* HV_chn(:) ); end
%     for parId1=1:numPar
%         for parId2=1:numPar
%             indc=1:numPar; indr=ones(1,numPar); indc( (indc==parId1)|(indc==parId2) )=[]; prob_val = HV_chn;
%             for ind3=1:numel(indc),prob_val = sum(prob_val,indc(ind3)); indr(indc(ind3))=siz(indc(ind3)); end
%             prob_val = repmat(prob_val,indr); prob_val = prob_val / sum(prob_val(:));
%             cov_mat(parId1,parId2) = sum((par(:,parId1)-mean_mat(parId1)) .* (par(:,parId2)-mean_mat(parId2)) .* prob_val(:)); 
%         end
%     end
%     % Standard Deviation
%     conf_values(channelId) = power(det(cov_mat),1/6);
%     max_prob_val(channelId) = sum(sum(sum(HV_chn(r_m-1:r_m+1,c_m-1:c_m+1,s_m-1:s_m+1))));
%   end
%   fileID = fopen(opts.entropysavename,'w'); display(max_prob_val);
%   fprintf(fileID,'ent: %3.3f\n' ,mean(max_prob_val));
%   fclose(fileID);
%   
  % compute the confidence of PEMS
  siz = size(HV); HV = HV + 1e-20;
  HV = HV ./ repmat(sum(sum(sum(HV,1),2),3),siz(1),siz(2),siz(3),1);
  HV = reshape(HV,prod(siz(1:3)),siz(4));
  entropies = -1*sum(HV.*log(HV),1);
  fileID = fopen(opts.entropysavename,'w');
  fprintf(fileID,'ent: %3.3f\n' ,mean(entropies));
  fclose(fileID);
  
end

%% Preprocessing function for the input images / The same procedure is
%%%% required for the training data
function funPreprocessing (op)

%%%% Resample to fixed isotropic resolution 
binaryresample  ='irtk_binaries/linux/resample';
system(sprintf('%s %s %s -linear -size %f %f %f',binaryresample,op.imagename,op.tmpimagename, op.pixel_spacing(1), op.pixel_spacing(2), op.pixel_spacing(3))); 

%%%% If it's the last stage - then perform the normalization based on the cropped image
if(op.stageId == 1)
    binaryregion = 'irtk_binaries/linux/region2';
    system(sprintf('%s %s %s -landmarks %s',binaryregion,op.tmpimagename,op.tmpimagename,op.vtkname)); 
end

%%%% Histogram Normalization ( not for CT images - other modalities )
binarynormalize ='irtk_binaries/linux/normalize';
if (op.stageId == 1), 
    referencename = op.referencecropname;
    system(sprintf('%s %s %s %s -piecewise',binarynormalize,referencename,op.tmpimagename,op.tmpimagename)); 
else
    referencename = op.referenceimagename;
    system(sprintf('%s %s %s %s -piecewise',binarynormalize,referencename,op.tmpimagename,op.tmpimagename)); 
end

%%%% Convert the images to NIFTI format (can be used to convert to NIFTI format)
binaryconvert   ='irtk_binaries/linux/convert';
system(sprintf('%s %s %s -ushort',binaryconvert,op.tmpimagename,op.tmpimagename)); 

end