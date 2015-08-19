__author__ = 'ozan'

import os
import subprocess

# Sbatch function to pass a command to slurm cluster
def sbatch(cmd, mem=5, c=10, n=1, queue='short', verbose=False, dryrun=False):
    sbatch_cmd = 'sbatch --mem='+str(mem)+'G' + ' -n ' + str(n) + ' -c '+str(c)+' -p '+queue+' --wrap="'+''.join(cmd)+'"'
    if verbose or dryrun:
        print sbatch_cmd
    if dryrun:
        return
    return subprocess.call(sbatch_cmd, shell=True)

def createFolder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

# Base Parameter List - Landmark Locations
slurm_ncores   = 8
slurm_nthreads = 1
slurm_memory   = 20
slurm_queue    = 'short'
source_dir     = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mri_edge_detect'
testdata_dir   = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mritestingdata'
modelname      = 'mriSecond_hier_exp_onlyClas'

input_img_dir  = testdata_dir + '/images'
ground_lm_dir  = testdata_dir + '/landmarks'
output_dir     = testdata_dir + '/pems_' + modelname; createFolder(output_dir)
results_dir    = testdata_dir + '/results';           createFolder(results_dir)

inputimages          = []
outputimages         = []
generatedLandmarks   = []
groundtruthLandmarks = []
for file in os.listdir(input_img_dir):
    if file.endswith('.nii.gz') and ('ed' in file):
        parsedfile = file.split('.nii.gz')
        inputimages.append(input_img_dir+'/'+file)
        outputimages.append(output_dir+'/'+file)
        groundtruthLandmarks.append(ground_lm_dir+'/'+parsedfile[0]+'.vtk')
        generatedLandmarks.append(output_dir+'/'+parsedfile[0]+'_lm.vtk')
inputimages          = sorted(inputimages)
outputimages         = sorted(outputimages)
groundtruthLandmarks = sorted(groundtruthLandmarks)
generatedLandmarks   = sorted(generatedLandmarks)

# Evaluation parameters
numSubjects     = min([1200,len(inputimages)])
distanceTxtFile = results_dir + '/distances_' + modelname + '.txt'
if os.path.exists(distanceTxtFile):
    os.remove(distanceTxtFile)

for ind in range(numSubjects):

  if not (os.path.exists(outputimages[ind])):

      # matlab multi-atlas script command
      cmd_pem   = '/usr/lib/matlab/R2014a/bin/matlab -nodesktop -nosplash -r \\"cd(\'{0}\'); addpath(\'{0}\'); edgesTestingDemo(\'{1}\',\'{2}\',\'{3}\'); quit;\\"'.format(source_dir,modelname,inputimages[ind],outputimages[ind])
      cmd_eval  = '/vol/medic02/users/oo2113/Build/irtk_master/bin/pevaluation {0} {1} -output {2}'.format( groundtruthLandmarks[ind], generatedLandmarks[ind], distanceTxtFile )
      cmd_slurm = cmd_pem + '; ' + cmd_eval
      # Run the proposed method (PEM) on SLURM (Multi-atlas segmentation)
      sbatch(cmd_slurm, mem=slurm_memory, n=slurm_nthreads, c=slurm_ncores, queue=slurm_queue, verbose=True, dryrun=False)
      #os.system(cmd_eval)

      '''
      cmd_pem   = '/usr/lib/matlab/R2014a/bin/matlab -nodesktop -nosplash -r \\"cd(\'\\\'\'{0}\'\\\'\'); edgesTestingDemo(\'\\\'\'{1}\'\\\'\',\'\\\'\'{2}\'\\\'\',\'\\\'\'{3}\'\\\'\'); quit;\\"'.format(source_dir,modelname,inputimages[ind],outputimages[ind])
      cmd_eval  = '/vol/medic02/users/oo2113/Build/irtk_master/bin/pevaluation {0} {1} -output {2}'.format( groundtruthLandmarks[ind], generatedLandmarks[ind], distanceTxtFile )
      cmd_pmk   = cmd_pem + '; ' + cmd_eval
      os.system('pmk.py --start --command "{0}"'.format(cmd_pmk))
      #os.system('{0}'.format(cmd_eval))
      '''
