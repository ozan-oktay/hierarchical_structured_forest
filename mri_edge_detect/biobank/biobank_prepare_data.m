function biobank_prepare_data()

%Add the libraries to MATLAB's path
restoredefaultpath;
setenv('LD_LIBRARY_PATH',[':/vol/medic02/users/oo2113/LOCAL/lib' ':/vol/medic02/users/oo2113/LOCAL/bin' ]);
setenv('PYTHONPATH',':/vol/bitbucket/oo2113/lib/python2.7/site-packages')

matlabpool(5);

%%%% Parameters
intensity_range_b = 0;
intensity_range_e = 1024;
voxelSize         = 1.25;
z_spacing         = 2.00;
label_type        = 'uchar';
image_type        = 'ushort';

%%%% Define the directories
addpath (genpath('/vol/medic02/users/oo2113/registration/src/matlab/denoise_cine'));
addpath (genpath('/homes/oo2113/workspace/matlab/nifti-read'));
addpath (genpath('/vol/biomedic/users/oo2113/str_hier_forest/toolbox'));

targetdir_training = '/vol/biomedic/users/oo2113/str_hier_forest_mri/biobankdata';
atlasdirectory     = '/vol/vipdata/data/biobank/cardiac/nifti';
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
mkdir(strcat(targetdir_training,'/landmarks'));
mkdir(strcat(targetdir_training,'/reference'));

%%%% Read the folder of atlas images 
d = dir(strcat(atlasdirectory,'/*')); 
atlasfoldernames = {d([d(:).isdir]).name}'; atlasfoldernames(ismember(atlasfoldernames,{'.','..'})) = [];

%%%% Define a reference image & rescale it and place it
system ((sprintf('%s %s %s',convertbin,ref,refimagename)));
system (sprintf('rescale %s %s "%d" "%d"',refimagename,refimagename,intensity_range_b,intensity_range_e));
copyfile(reflm,reflandmarkname);

%%%% Pick an atlas in a for loop 
parfor patientId = 1 : numel (atlasfoldernames)
  
  %%%% Check if the images already exist 
  patientName = atlasfoldernames{patientId};
  img_exist   = isempty(dir(strcat(targetdir_training,'/images',sprintf('/%s_*',patientName))));
  if (~img_exist)
    continue;
  end
  fprintf ('processing images for patient %d\n',patientId);
  
  %%%% Create the temporary folder (specific to workerId)
  t = getCurrentTask(); 
  tmpfoldername = sprintf('tmp%d',t.ID);
  targetdir_training_tmp = strcat(targetdir_training,'/',tmpfoldername);
  mkdir(targetdir_training_tmp);
  
  %%%% Remove first - Convert it to .nii.gz uint16
  atlasfolderdir = strcat(atlasdirectory,'/',atlasfoldernames{patientId});
  volume4dname   = strcat(targetdir_training_tmp,'/cine_sax_4D.nii.gz');

  %%%% Determine the ED and ES frames
  ed_frame = 0; system(horzcat(convertbin,' ',atlasfolderdir,'/cine_sax.nii.gz',' ',volume4dname,' -ushort'));
      
  %%%% Rescale the volume image 
  system(sprintf('rescale_real %s %s "%d" "%d"',volume4dname,volume4dname,intensity_range_b,intensity_range_e));
  
  %%%% Normalize the volume image
  system(sprintf('normalize %s %s %s -Tp 0 -Sp 0 -piecewise',refimagename,volume4dname,volume4dname));
      
  %%%% Clear & Split the 3D+T into separate frames
  system(horzcat('rm -rf ',targetdir_training_tmp,'/frame*.nii.gz'));
  system(horzcat('splitvolume',' ',volume4dname,' ',targetdir_training_tmp,'/frame',' ','-sequence'));
  tmp_edimgname   = strcat(targetdir_training_tmp,sprintf('/frame%2.2d.nii.gz',ed_frame));
  
  %%%% Read Voxel spacing and resize inplane to same voxel level (both image and label)
  system(sprintf('resample %s %s -size %g %g %g -linear -nne',tmp_edimgname,tmp_edimgname,voxelSize,voxelSize,z_spacing));
  target_imgnameED = strcat(targetdir_training,'/images',sprintf('/%s_ed%d_3D.nii.gz',patientName,ed_frame));
  system (horzcat(convertbin,' ',tmp_edimgname,' ',target_imgnameED,' ','-',image_type));

  %%%% Save the VTK Landmarks as a vtk file 
  target_landmarksnameED_vtk = strcat(targetdir_training,'/landmarks',sprintf('/%s_ed%d_3D.vtk',patientName,ed_frame));
  target_landmarksnameED_txt = strcat(targetdir_training,'/landmarks',sprintf('/%s_ed%d_3D.txt',patientName,ed_frame));
  copyfile(strcat(atlasfolderdir,'/landmarks.vtk'),target_landmarksnameED_vtk);
  
  %%%% Reset the headers - After saving the landmarks in image coordinate space
  system(horzcat('headertool',' ',target_imgnameED,' ',tmp_edimgname,' ','-reset'));
  system(horzcat('ptransformation',' ',target_landmarksnameED_vtk,' ',target_landmarksnameED_vtk,' -target ',target_imgnameED,' -source ',tmp_edimgname));
  system(horzcat('headertool',' ',tmp_edimgname,' ',target_imgnameED,' ','-orientation',' ','1 0 0 0 -1 0 0 0 -1'));
  system(horzcat('ptransformation',' ',target_landmarksnameED_vtk,' ',target_landmarksnameED_vtk,' -target ',tmp_edimgname,' -source ',target_imgnameED));
  system(horzcat('transformation',' ',target_imgnameED,' ',target_imgnameED,' ','-target',' ',tmp_edimgname,' ','-linear -matchInputType'));

  %%%% Convert the VTK Landmarks to a text file
  system(horzcat('vtk2txt.py',' ','--input',' ',target_landmarksnameED_vtk,' ','--output',' ',target_landmarksnameED_txt,' ','--imagecoord true',' ','--targetimage',' ',target_imgnameED));

  %%%% Clean the temporary directory
  system(horzcat('rm -rf ',targetdir_training_tmp));
  
end

system(horzcat('/vol/biomedic/users/oo2113/str_hier_forest_mri/mri_edge_detect/python/computeOrientations.py --indir ',sprintf('%s',targetdir_training)));
matlabpool close;

