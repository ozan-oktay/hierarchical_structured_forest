__author__ = 'ozan'

import os
import subprocess

def mycreatedir (l_directory):
    if not os.path.exists(l_directory):
        os.makedirs(l_directory)

def sbatch(cmd, mem=5, c=10, n=1, queue='short', verbose=False, dryrun=False):
    sbatch_cmd = 'sbatch --mem='+str(mem)+'G' + ' -n ' + str(n) + ' -c '+str(c)+' -p '+queue+' --wrap="'+''.join(cmd)+'"'

    if verbose or dryrun:
        print sbatch_cmd

    if dryrun:
        return

    return subprocess.call(sbatch_cmd, shell=True)

slurm_ncores   = 8
slurm_nthreads = 1
slurm_memory   = 100
slurm_queue    = 'interactive'
source_dir     = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mri_edge_detect'

# matlab multi-atlas script command
cmd_pem = '/usr/lib/matlab/R2015a/bin/matlab -nodesktop -nosplash -r \\"addpath(\'{0}\'); main_prepare_data_large(); quit;\\"'.format(source_dir)
# Run the proposed method (PEM) on SLURM (Multi-atlas segmentation)
sbatch(cmd_pem, mem=slurm_memory, n=slurm_nthreads, c=slurm_ncores, queue=slurm_queue, verbose=True, dryrun=False)
