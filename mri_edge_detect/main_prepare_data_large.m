function main_prepare_data_large()

%Add the libraries to MATLAB's path
restoredefaultpath;
setenv('LD_LIBRARY_PATH',[':/vol/medic02/users/oo2113/LOCAL/lib' ':/vol/medic02/users/oo2113/LOCAL/bin' ]);
setenv('PYTHONPATH',':/vol/bitbucket/oo2113/lib/python2.7/site-packages')

pool=parpool(4); 
 
%%%% Parameters
intensity_range_b = 0;
intensity_range_e = 1024;
voxelSize         = 1.25;
z_spacing         = 2.00;
dofrange          = [-75,25];
label_type        = 'uchar';
image_type        = 'ushort';
selectedLabels    = [1,2,4];

%%%% Define the directories
setenv  ('LD_LIBRARY_PATH', '');
addpath (genpath('/vol/medic02/users/oo2113/registration/src/matlab/denoise_cine'));
addpath (genpath('/homes/oo2113/workspace/matlab/nifti-read'));
addpath (genpath('/vol/biomedic/users/oo2113/str_hier_forest_mri/toolbox'));

targetdir_training = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mritrainingdata_sec_large2';
atlasdirectory     = '/vol/medic02/users/wbai/data/3d_data_1000';
reffolder          = '/vol/vipdata/data/cardiac/3datlas/1_FINAL/';

ref                = strcat(reffolder,'/lvsa_ED.gipl');
reflm              = strcat(reffolder,'/landmarks.vtk');
refimagename       = strcat(targetdir_training,'/reference/ref_image.nii.gz');
reflandmarkname    = strcat(targetdir_training,'/reference/ref_lm.vtk');
convertbin         = '/vol/medic02/users/oo2113/Build/irtk/bin/convert';

%%%% Create training data folders 
mkdir(targetdir_training);
mkdir(strcat(targetdir_training,'/tmp'));
mkdir(strcat(targetdir_training,'/images'));
mkdir(strcat(targetdir_training,'/groundtruth'));
mkdir(strcat(targetdir_training,'/landmarks'));
mkdir(strcat(targetdir_training,'/reference'));
mkdir(strcat(targetdir_training,'/boundary'));

%%%% Read the folder of atlas images 
d         = [dir(strcat(atlasdirectory,'/14*'))];
d_first   = [dir(strcat(atlasdirectory,'/14A*'));dir(strcat(atlasdirectory,'/14B*'));dir(strcat(atlasdirectory,'/14C*'));];
d_testing = [dir(strcat(atlasdirectory,'/14D*'));dir(strcat(atlasdirectory,'/14E*'));dir(strcat(atlasdirectory,'/14F*'));dir(strcat(atlasdirectory,'/14G*'))];
d_excepts = [d_first;d_testing]; 

atlasfoldernames    = {d([d(:).isdir]).name}';                 atlasfoldernames(ismember(atlasfoldernames,{'.','..'}))       = [];
unwantedfoldernames = {d_excepts([d_excepts(:).isdir]).name}'; unwantedfoldernames(ismember(unwantedfoldernames,{'.','..'})) = [];
atlasfoldernames    = setdiff(atlasfoldernames,unwantedfoldernames);

%%%% Define a reference image & rescale it and place it
system ((sprintf('%s %s %s',convertbin,ref,refimagename)));
system (sprintf('rescale %s %s "%d" "%d"',refimagename,refimagename,intensity_range_b,intensity_range_e));
copyfile(reflm,reflandmarkname);

%%%% Pick an atlas in a for loop 
parfor patientId = 1 : numel (atlasfoldernames)
  
  %%%% Check if the images already exist 
  patientName    = atlasfoldernames{patientId};
  patientNameCon = strcat(atlasfoldernames{patientId},'_con');
  img_exist      = isempty(dir(strcat(targetdir_training,'/images',sprintf('/%s_*',patientName))));
  gt_exist       = isempty(dir(strcat(targetdir_training,'/groundtruth',sprintf('/%s_*',patientName))));
  boun_exist     = isempty(dir(strcat(targetdir_training,'/boundary',sprintf('/%s_*',patientName))));
  if (~img_exist && ~gt_exist && ~boun_exist)
    continue;
  end
  fprintf ('processing images for patient %d - %s\n',patientId,patientName);
  
  %%%% Create the temporary folder (specific to workerId)
  t = getCurrentTask(); 
  tmpfoldername = sprintf('tmp%d',t.ID);
  targetdir_training_tmp = strcat(targetdir_training,'/',tmpfoldername);
  mkdir(targetdir_training_tmp);
  
  %%%% Remove first - Convert it to .nii.gz uint16
  atlasfolderdir = strcat(atlasdirectory,'/',atlasfoldernames{patientId});
  volume4dname   = strcat(targetdir_training_tmp,'/volume4d.nii.gz');
  edimgname      = strcat(targetdir_training_tmp,'/ed.nii.gz');
  esimgname      = strcat(targetdir_training_tmp,'/es.nii.gz');
  edlabelname    = strcat(targetdir_training_tmp,'/edseg.nii.gz');
  eslabelname    = strcat(targetdir_training_tmp,'/esseg.nii.gz');

  %%%% Determine the ED and ES frames
  ed_frame = -1; es_frame = -1;
  system(horzcat(convertbin,' ',atlasfolderdir,'/lvsa.nii.gz',' ',volume4dname,' -ushort'));
  system(horzcat(convertbin,' ',atlasfolderdir,'/lvsa_ED.nii.gz',' ',edimgname,' -ushort'));
  system(horzcat(convertbin,' ',atlasfolderdir,'/lvsa_ES.nii.gz',' ',esimgname,' -ushort'));
  system(horzcat(convertbin,' ',atlasfolderdir,'/segmentation_ED.nii.gz',' ',edlabelname,' -uchar'));
  system(horzcat(convertbin,' ',atlasfolderdir,'/segmentation_ES.nii.gz',' ',eslabelname,' -uchar'));
  vol4d_img = load_untouch_nii(volume4dname); vol4d_img = vol4d_img.img;
  ed_img    = load_untouch_nii(edimgname);       ed_img    = ed_img.img;
  es_img    = load_untouch_nii(esimgname);       es_img    = es_img.img;
  
  if (~isequal(size(ed_img),size(vol4d_img(:,:,:,1))) || ~isequal(size(es_img),size(vol4d_img(:,:,:,1))))
    continue;
  end
  
  for frameId = 1 : size(vol4d_img,4)
    diff_ed = (vol4d_img (:,:,:,frameId) - ed_img);
    diff_es = (vol4d_img (:,:,:,frameId) - es_img);
    if (diff_ed==0), ed_frame = frameId -1; end % C convention start from 0...n
    if (diff_es==0), es_frame = frameId -1; end % C convention start from 0...n
  end
  if ((ed_frame < 0) || (es_frame < 0)), error('couldnt find the corresponding ed or es frame'); end;
      
  %%%% Rescale the volume image 
  system(sprintf('rescale_real %s %s "%d" "%d"',volume4dname,volume4dname,intensity_range_b,intensity_range_e));
  
  %%%% Normalize the volume image
  system(sprintf('normalize %s %s %s -Tp 0 -Sp 0 -piecewise',refimagename,volume4dname,volume4dname));

  %%%% Read Voxel spacing and resize inplane to same voxel level (both image and label)
  system(sprintf('resample %s %s -size %g %g %g -linear -nne',volume4dname,volume4dname,voxelSize,voxelSize,z_spacing));
  
  % Map the labels into the same space of 4D Image 
  system(sprintf('transformation %s %s -dofin "Id" -target %s -nn -matchInputType',edlabelname,edlabelname,volume4dname));
  
  % Selected the desired label IDs (e.g. endocardium only)
  system(sprintf('selectLabelId.py --inputname %s --outputname %s -l %s',edlabelname,edlabelname,sprintf('%d ',selectedLabels)));
  
  %%%% Clear & Split the 3D+T into separate frames
  system(horzcat('rm -rf ',targetdir_training_tmp,'/frame*.nii.gz'));
  system(horzcat('splitvolume',' ',volume4dname,' ',targetdir_training_tmp,'/frame',' ','-sequence'));
  tmp_edimgname = strcat(targetdir_training_tmp,sprintf('/frame%2.2d.nii.gz',ed_frame));
  
  target_imgnameED   = strcat(targetdir_training,'/images',sprintf('/%s_ed%d_3D.nii.gz',patientName,ed_frame));
  target_labelnameED = strcat(targetdir_training,'/groundtruth',sprintf('/%s_ed%d_3D.nii.gz',patientName,ed_frame));
  
  system (horzcat(convertbin,' ',tmp_edimgname,' ',target_imgnameED,' ','-',image_type));
  system (horzcat(convertbin,' ',edlabelname,' '  ,target_labelnameED,' ','-',label_type));

  %%%% Generate Boundaries 
  target_boundnameED    = strcat(targetdir_training,'/boundary',sprintf('/%s_ed%d_3D.nii.gz',patientName,ed_frame));
  system(horzcat('sobelEdgeDetector.py',' ','--img',' ',target_labelnameED,' ','--output',' ',target_boundnameED));
    
  %%%% Save the VTK Landmarks as text file 
  target_landmarksnameED_txt = strcat(targetdir_training,'/landmarks',sprintf('/%s_ed%d_3D.txt',patientName,ed_frame));
  target_landmarksnameED_vtk = strcat(targetdir_training,'/landmarks',sprintf('/%s_ed%d_3D.vtk',patientName,ed_frame));
  system(horzcat('vtk2txt.py',' ','--input',' ',atlasfolderdir,'/landmarks.vtk',' ','--output',' ',target_landmarksnameED_txt,' ','--imagecoord true',' ','--targetimage',' ',target_imgnameED));
  copyfile(strcat(atlasfolderdir,'/landmarks.vtk'),target_landmarksnameED_vtk);
  
  %%%% Reset the headers - After saving the landmarks in image coordinate space
  system(horzcat('headertool',' ',target_imgnameED,' ',target_imgnameED,' ','-reset'));
  system(horzcat('headertool',' ',target_labelnameED,' ',target_labelnameED,' ','-reset'));
  
  %%%% Rotate the images, the label maps and generate the boundary images
  target_imgnameEDCon           = strcat(targetdir_training,'/images',sprintf('/%s_ed%d_3D.nii.gz',patientNameCon,ed_frame));
  target_labelnameEDCon         = strcat(targetdir_training,'/groundtruth',sprintf('/%s_ed%d_3D.nii.gz',patientNameCon,ed_frame));
  target_boundnameEDCon         = strcat(targetdir_training,'/boundary',sprintf('/%s_ed%d_3D.nii.gz',patientNameCon,ed_frame));
  target_landmarksnameEDCon_vtk = strcat(targetdir_training,'/landmarks',sprintf('/%s_ed%d_3D.vtk',patientNameCon,ed_frame));
  target_landmarksnameEDCon_txt = strcat(targetdir_training,'/landmarks',sprintf('/%s_ed%d_3D.txt',patientNameCon,ed_frame));
  randomDofname                 = strcat(tempname,'.dof.gz');
  system(sprintf('dofcreate %s -rz %f',randomDofname,randi(dofrange,1,1)));
  system(sprintf('transformation %s %s -dofin %s -target %s -linear',target_imgnameED,target_imgnameEDCon,randomDofname,target_imgnameED));
  system(sprintf('transformation %s %s -dofin %s -target %s -nn ',   target_labelnameED,target_labelnameEDCon,randomDofname,target_labelnameED));
  system(horzcat('sobelEdgeDetector.py',' ','--img',' ',target_labelnameEDCon,' ','--output',' ',target_boundnameEDCon));
  system(sprintf('ptransformation %s %s -invert -dofin %s',target_landmarksnameED_vtk,target_landmarksnameEDCon_vtk,randomDofname))
  system(horzcat('vtk2txt.py',' ','--input',' ',target_landmarksnameEDCon_vtk,' ','--output',' ',target_landmarksnameEDCon_txt,' ','--imagecoord true',' ','--targetimage',' ',target_imgnameEDCon));

  %%%% Clean the temporary directory
  system(horzcat('rm -rf ',targetdir_training_tmp));
  system(horzcat('rm -rf ',randomDofname));
  
end

system(horzcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/mri_edge_detect/python/computeOrientations.py --indir ',sprintf('%s',targetdir_training)));
poolobj = gcp('nocreate'); delete(poolobj); 
