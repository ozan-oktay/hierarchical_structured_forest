__author__ = 'ozan'

import os
import subprocess

# Sbatch function to pass a command to slurm cluster
def createFolder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

# Base Parameter List - Landmark Locations
source_dir     = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mri_edge_detect'
testdata_dir   = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mribiobankdata'

ground_lm_dir1  = testdata_dir + '/landmarks'
ground_lm_dir2  = testdata_dir + '/landmarks_ozan'
results_dir     = testdata_dir + '/results';           createFolder(results_dir)

groundtruthLandmarks1 = []
groundtruthLandmarks2 = []

for file in os.listdir(ground_lm_dir2):
    if file.endswith('.vtk') and ('ed' in file):
        parsedfile = file.split('.vtk')
        groundtruthLandmarks1.append(ground_lm_dir1+'/'+parsedfile[0]+'.vtk')
        groundtruthLandmarks2.append(ground_lm_dir2+'/'+parsedfile[0]+'.vtk')

groundtruthLandmarks1 = sorted(groundtruthLandmarks1)
groundtruthLandmarks2 = sorted(groundtruthLandmarks2)
numSubjects           = len(groundtruthLandmarks2)

# Evaluation parameters
distanceTxtFile  = results_dir + '/distances_interUser.txt'
if os.path.exists(distanceTxtFile):  os.remove(distanceTxtFile)

for ind in range(numSubjects):

    cmd_eval  = '/vol/medic02/users/oo2113/Build/irtk_master/bin/pevaluation {0} {1} -output {2}'.format( groundtruthLandmarks1[ind], groundtruthLandmarks2[ind], distanceTxtFile )
    os.system(cmd_eval)
