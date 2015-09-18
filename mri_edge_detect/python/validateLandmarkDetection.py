__author__ = 'ozan'

import os
import subprocess

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

def evalPose(groundtruthPoseFile,generatedPoseFile,poseErrorTxt):
  import irtk
  import re
  # Read the ground-truth file
  grndDof   = irtk.AffineTransformation(filename=groundtruthPoseFile)
  grndArray = [grndDof.rz,grndDof.rx,grndDof.ry]

  # Read the generated pose text file
  txtfile     = open(generatedPoseFile, 'r')
  txtline     = txtfile.readline(); computedMean = re.findall(r"[-+]?\d*\.\d+|\d+", txtline)
  txtline     = txtfile.readline(); computedStd  = re.findall(r"[-+]?\d*\.\d+|\d+", txtline)
  txtfile.close()

  # Write the values
  with open(poseErrorTxt, "a") as myfile:
    myfile.write('file: {0} grnd: {1} cmpt: {2} conf(std): {3}\n'.format(inputimages[ind],grndArray[0],computedMean[0],computedStd[0]))
    myfile.close()

# Base Parameter List - Landmark Locations
slurm_ncores   = 8
slurm_nthreads = 1
slurm_memory   = 25
slurm_queue    = 'short'
source_dir     = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mri_edge_detect'
testdata_dir   = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mribiobankdata'
modelname      = 'mriSecond_hier_Mix3'

input_img_dir  = testdata_dir + '/images'
ground_lm_dir  = testdata_dir + '/landmarks'
ground_dof_dir = testdata_dir + '/dofs'
output_dir     = testdata_dir + '/pems_' + modelname; createFolder(output_dir)
results_dir    = testdata_dir + '/results';           createFolder(results_dir)

inputimages          = []
outputimages         = []
generatedLandmarks   = []
groundtruthLandmarks = []
generatedPoseFiles   = []
groundtruthPoseFiles = []

for file in os.listdir(input_img_dir):
    if file.endswith('.nii.gz') and ('ed' in file):
        parsedfile = file.split('.nii.gz')
        inputimages.append(input_img_dir+'/'+file)
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

# Evaluation parameters
numSubjects      = min([1200,len(inputimages)])
print ('number of subjects is {0}'.format(numSubjects))
distanceTxtFile  = results_dir + '/distances_'       + modelname + '.txt'
poseErrorTxtFile = results_dir + '/poseEstimations_' + modelname + '.txt'
if os.path.exists(distanceTxtFile):  os.remove(distanceTxtFile)
if os.path.exists(poseErrorTxtFile): os.remove(poseErrorTxtFile)

for ind in range(numSubjects):

  if os.path.exists(outputimages[ind]):

    cmd_eval  = '/vol/medic02/users/oo2113/Build/irtk_master/bin/pevaluation {0} {1} -output {2}'.format( groundtruthLandmarks[ind], generatedLandmarks[ind], distanceTxtFile )
    os.system(cmd_eval)
    #if os.path.exists(generatedPoseFiles[ind]): evalPose(groundtruthPoseFiles[ind],generatedPoseFiles[ind],poseErrorTxtFile)

  else:

      # matlab multi-atlas script command
      cmd_pem   = '/usr/lib/matlab/R2015a/bin/matlab -nodesktop -nosplash -r \\"cd(\'{0}\'); addpath(\'{0}\'); edgesTestingDemo(\'{1}\',\'{2}\',\'{3}\'); quit;\\"'.format(source_dir,modelname,inputimages[ind],output_dir)
      cmd_eval  = '/vol/medic02/users/oo2113/Build/irtk_master/bin/pevaluation {0} {1} -output {2}'.format( groundtruthLandmarks[ind], generatedLandmarks[ind], distanceTxtFile )
      logname   = '/vol/bitbucket/oo2113/tmp/logfile.out'
      cmd_slurm = cmd_pem + '; ' + cmd_eval
      # Run the proposed method (PEM) on SLURM (Multi-atlas segmentation)
      sbatch(cmd_slurm, logname, mem=slurm_memory, n=slurm_nthreads, c=slurm_ncores, queue=slurm_queue, verbose=True, dryrun=False)

      '''
      cmd_pem   = '/usr/lib/matlab/R2015a/bin/matlab -nodesktop -nosplash -r \\"cd(\'\\\'\'{0}\'\\\'\'); edgesTestingDemo(\'\\\'\'{1}\'\\\'\',\'\\\'\'{2}\'\\\'\',\'\\\'\'{3}\'\\\'\'); quit;\\"'.format(source_dir,modelname,inputimages[ind],outputimages[ind])
      cmd_eval  = '/vol/medic02/users/oo2113/Build/irtk_master/bin/pevaluation {0} {1} -output {2}'.format( groundtruthLandmarks[ind], generatedLandmarks[ind], distanceTxtFile )
      cmd_pmk   = cmd_pem + '; ' + cmd_eval
      os.system('pmk.py --start --command "{0}"'.format(cmd_pmk))
      #os.system('{0}'.format(cmd_eval))
      '''

