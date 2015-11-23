__author__ = 'ozan'

import os
import subprocess
import numpy as np
import tempfile as tp

# Sbatch function to pass a command to slurm cluster
def sbatch(cmd, logname, mem=5, c=10, n=1, queue='short', verbose=False, dryrun=False):
    sbatch_cmd = 'sbatch --mem='+str(mem)+'G' + ' -n ' + str(n) + ' -o ' +str(logname)+ ' -c '+str(c)+' -p '+queue+' --wrap="'+''.join(cmd)+'"'
    if verbose or dryrun:
        print sbatch_cmd
    if dryrun:
        return
    return subprocess.call(sbatch_cmd, shell=True)

def createFolder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def evalPose(imagename,groundtruthPoseFile,generatedPoseFile,poseErrorTxt):
  import irtk
  import re

  # Read the ground-truth file
  grndDof   = irtk.AffineTransformation(filename=groundtruthPoseFile)
  affine_trans = np.array(grndDof.matrix())
  determinant  = np.linalg.det(affine_trans)
  grndScale    = np.power(determinant,1.0/3.0)
  grndRot      = grndDof.rz

  # Read the generated pose text file
  txtfile     = open(generatedPoseFile, 'r')
  txtline     = txtfile.readline(); computedMean = re.findall(r"[-+]?\d*\.\d+|\d+", txtline)
  txtline     = txtfile.readline(); computedStd  = re.findall(r"[-+]?\d*\.\d+|\d+", txtline)
  txtfile.close()

  # Write the values
  with open(poseErrorTxt, "a") as myfile:
    myfile.write('file: {0} '.format(imagename))
    myfile.write('grnd: {0} {1} '.format(grndRot,grndScale))
    myfile.write('cmpt: ')
    for mean in computedMean:
      myfile.write('{0} '.format(mean))
    myfile.write('conf(std): ')
    for std in computedStd:
      myfile.write('{0} '.format(std))
    myfile.write('\n')
    myfile.close()

# Base Parameter List - Landmark Locations
slurm_ncores   = 8
slurm_nthreads = 1
slurm_memory   = 25
slurm_queue    = 'short'
slurm_logname  = '/vol/bitbucket/oo2113/tmp/logfile.out'
base_dir       = '/vol/medic02/users/oo2113/str_hier_forest_mri'
source_dir     = base_dir + '/mri_edge_detect'
testdata_dir   = base_dir + '/mritestingdata'
modelname      = 'new_prn016_shape_large'
n_atlases      = 3

input_img_dir    = testdata_dir + '/images'
ground_lm_dir    = testdata_dir + '/landmarks'
ground_label_dir = testdata_dir + '/groundtruth'
ground_dof_dir   = testdata_dir + '/dofs'
output_dir       = testdata_dir + '/pems_' + modelname; createFolder(output_dir)
results_dir      = testdata_dir + '/results';           createFolder(results_dir)
reference_dir    = testdata_dir + '/reference'
atlas_dir        = base_dir + '/segmentation/atlases'

inputimages          = []
outputimages         = []
generatedLandmarks   = []
groundtruthLandmarks = []
generatedPoseFiles   = []
groundtruthPoseFiles = []
groundtruthLabels    = []

for file in os.listdir(input_img_dir):
    if file.endswith('.nii.gz') and ('ed' in file):
        parsedfile = file.split('.nii.gz')
        inputimages.append(input_img_dir+'/'+file)
        groundtruthLabels.append(ground_label_dir+'/'+file)
        outputimages.append(output_dir+'/'+parsedfile[0]+'_pem.nii.gz')
        groundtruthLandmarks.append(ground_lm_dir+'/'+parsedfile[0]+'.vtk')
        generatedLandmarks.append(output_dir+'/'+parsedfile[0]+'_lm.vtk')
        groundtruthPoseFiles.append(ground_dof_dir+'/'+parsedfile[0]+'.dof.gz')
        generatedPoseFiles.append(output_dir+'/'+parsedfile[0]+'_pose.txt')

inputimages          = sorted(inputimages)
outputimages         = sorted(outputimages)
groundtruthLandmarks = sorted(groundtruthLandmarks)
generatedLandmarks   = sorted(generatedLandmarks)
groundtruthPoseFiles = sorted(groundtruthPoseFiles)
generatedPoseFiles   = sorted(generatedPoseFiles)
groundtruthLabels    = sorted(groundtruthLabels)

# Evaluation parameters
numSubjects      = min([500000,len(inputimages)])
print ('number of subjects is {0}'.format(numSubjects))
distanceTxtFile      = results_dir + '/distances_'         + modelname + '.txt'
poseErrorTxtFile     = results_dir + '/poseEstimations_'   + modelname + '.txt'
rigidRegErrorTxtFile = results_dir + '/rigidRegistration_' + modelname + '.txt'
if os.path.exists(distanceTxtFile):  os.remove(distanceTxtFile)
if os.path.exists(poseErrorTxtFile): os.remove(poseErrorTxtFile)
if os.path.exists(rigidRegErrorTxtFile): os.remove(rigidRegErrorTxtFile)

for ind in range(numSubjects):

  if os.path.exists(outputimages[ind]):

    # landmark error evaluation
    cmd_eval  = '/vol/medic02/users/oo2113/Build/irtk_master/bin/pevaluation {0} {1} -output {2};'.format( groundtruthLandmarks[ind], generatedLandmarks[ind], distanceTxtFile )
    # rigid-registration error evaluation
    for ind2 in range(1,n_atlases+1):
        tempdof1  = tp.NamedTemporaryFile(delete = True,suffix='.dof.gz'); cmd_eval += ' /vol/medic02/users/oo2113/Build/irtk/bin/prreg {0} {1} -dofout {2};'.format(atlas_dir+'/Atlas{0:02d}/landmarks.vtk'.format(ind2),generatedLandmarks[ind],tempdof1.name)
        tempdof2  = tp.NamedTemporaryFile(delete = True,suffix='.dof.gz'); cmd_eval += ' /vol/medic02/users/oo2113/Build/irtk/bin/prreg {0} {1} -dofout {2};'.format(atlas_dir+'/Atlas{0:02d}/landmarks.vtk'.format(ind2),groundtruthLandmarks[ind],tempdof2.name)
        cmd_eval  += ' /vol/medic02/users/oo2113/Build/irtk_master/bin/revaluation {0} {1} {2} -output {3} -vtk {4};'.format(groundtruthLabels[ind],tempdof1.name,tempdof2.name,rigidRegErrorTxtFile,groundtruthLandmarks[ind])
    os.system(cmd_eval)
    if os.path.exists(generatedPoseFiles[ind]): evalPose(groundtruthLabels[ind],groundtruthPoseFiles[ind],generatedPoseFiles[ind],poseErrorTxtFile)

  else:

      # matlab multi-atlas script command
      cmd_pem   = '/usr/lib/matlab/R2015a/bin/matlab -nodesktop -nosplash -r \\"cd(\'{0}\'); addpath(\'{0}\'); edgesTestingDemo(\'{1}\',\'{2}\',\'{3}\'); quit;\\"'.format(source_dir,modelname,inputimages[ind],output_dir)

      # landmark error evaluation
      cmd_eval  = '/vol/medic02/users/oo2113/Build/irtk_master/bin/pevaluation {0} {1} -output {2};'.format( groundtruthLandmarks[ind], generatedLandmarks[ind], distanceTxtFile )

      # rigid-registration error evaluation
      for ind2 in range(1,n_atlases+1):
        tempdof1=tp.NamedTemporaryFile(delete = True,suffix='.dof.gz'); cmd_eval += ' /vol/medic02/users/oo2113/Build/irtk/bin/prreg {0} {1} -dofout {2};'.format(atlas_dir+'/Atlas{0:02d}/landmarks.vtk'.format(ind2),generatedLandmarks[ind],tempdof1.name)
        tempdof2=tp.NamedTemporaryFile(delete = True,suffix='.dof.gz'); cmd_eval += ' /vol/medic02/users/oo2113/Build/irtk/bin/prreg {0} {1} -dofout {2};'.format(atlas_dir+'/Atlas{0:02d}/landmarks.vtk'.format(ind2),groundtruthLandmarks[ind],tempdof2.name)
        cmd_eval += ' /vol/medic02/users/oo2113/Build/irtk_master/bin/revaluation {0} {1} {2} -output {3} -vtk {4};'.format(groundtruthLabels[ind],tempdof1.name,tempdof2.name,rigidRegErrorTxtFile,groundtruthLandmarks[ind])

      # Run the proposed method (PEM) on SLURM (Multi-atlas segmentation)
      cmd_slurm = cmd_pem + '; ' + cmd_eval
      sbatch(cmd_slurm, slurm_logname, mem=slurm_memory, n=slurm_nthreads, c=slurm_ncores, queue=slurm_queue, verbose=True, dryrun=False)


      '''
      cmd_pem   = '/usr/lib/matlab/R2015a/bin/matlab -nodesktop -nosplash -r \\"cd(\'\\\'\'{0}\'\\\'\'); edgesTestingDemo(\'\\\'\'{1}\'\\\'\',\'\\\'\'{2}\'\\\'\',\'\\\'\'{3}\'\\\'\'); quit;\\"'.format(source_dir,modelname,inputimages[ind],outputimages[ind])
      cmd_eval  = '/vol/medic02/users/oo2113/Build/irtk_master/bin/pevaluation {0} {1} -output {2}'.format( groundtruthLandmarks[ind], generatedLandmarks[ind], distanceTxtFile )
      cmd_pmk = cmd_pem + '; ' + cmd_eval + '; ' + cmd_eval2 + '; ' + cmd_eval3 + '; ' + cmd_eval4
      os.system('pmk.py --start --command "{0}"'.format(cmd_pmk))
      #os.system('{0}'.format(cmd_eval))
      '''

