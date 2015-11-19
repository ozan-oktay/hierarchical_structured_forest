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
slurm_queue = 'long'
source_dir = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mri_edge_detect/'
modelname=['prn016']
pr=[[0.16,0.42,0.42]]

for i in range(len(modelname)):
  # matlab multi-atlas script command
  cmd_pem = 'matlab -nodesktop -nosplash -r \\"addpath(\'{0}\'); edgesTrainingDemo(); quit;\\"'.format(source_dir,modelname[i],pr[i][0],pr[i][1],pr[i][2])
  
  # Run the proposed method (PEM) on SLURM (Multi-atlas segmentation)
  sbatch(cmd_pem, mem=slurm_memory, n=slurm_nthreads, c=slurm_ncores, queue=slurm_queue, verbose=True, dryrun=False)
