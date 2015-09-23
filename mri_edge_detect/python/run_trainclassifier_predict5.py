__author__ = 'ozan'

import os
import subprocess


def sbatch(cmd, mem=5, c=10, n=1, queue='short', verbose=False, dryrun=False):
    sbatch_cmd = 'sbatch --mem='+str(mem)+'G' + ' -n ' + str(n) + ' -c '+str(c)+' -p '+queue+' --wrap="'+''.join(cmd)+'"'

    if verbose or dryrun:
        print sbatch_cmd

    if dryrun:
        return

    return subprocess.call(sbatch_cmd, shell=True)

slurm_ncores = 1
slurm_nthreads = 8
slurm_memory = 250
slurm_queue = 'interactive'
source_dir = '/vol/medic02/users/oo2113/str_hier_forest_mri/mri_edge_detect/'
nodeProb=[[0.06,0.47,0.47],[0.1,0.45,0.45],[0.2,0.4,0.4],[0.333,0.333,0.333]];
modelname=['mymodel006','mymodel010','mymodel020','mymodel033']

for i in range(len(modelname)):
  # matlab multi-atlas script command
  cmd_pem = 'matlab -nodesktop -nosplash -r \\"addpath(\'{0}\'); edgesTrainingDemo([{1},{2},{3}],\'{4}\'); quit;\\"'.format(source_dir,nodeProb[i][0],nodeProb[i][1],nodeProb[i][2],modelname[i])
  
  # Run the proposed method (PEM) on SLURM (Multi-atlas segmentation)
  sbatch(cmd_pem, mem=slurm_memory, n=slurm_nthreads, c=slurm_ncores, queue=slurm_queue, verbose=True, dryrun=False)