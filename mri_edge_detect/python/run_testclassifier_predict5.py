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
slurm_memory   = 20
slurm_queue    = 'short'
source_dir     = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mri_edge_detect'
input_dir      = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mritestingdata/images'
output_dir     = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mritestingdata/pems_reg'; mycreatedir(output_dir)
modelname      = 'mri3Dmodel_reg_t50'

inputimages=[]
outputimages=[]

for file in os.listdir(input_dir):
    if file.endswith('.nii.gz') and ('ed' in file):
        parsedfile = file.split('.nii.gz')
        inputimages.append(input_dir+'/'+file)
        outputimages.append(output_dir+'/'+parsedfile[0]+'_pem.nii.gz')

for ind in range(len(inputimages)):

  # matlab multi-atlas script command
  #cmd_pem = '/usr/lib/matlab/R2014a/bin/matlab -nodesktop -nosplash -r \\"addpath(\'{0}\'); edgesTestingDemo(\'{1}\',\'{2}\',\'{3}\'); quit;\\"'.format(source_dir,modelname,inputimages[ind],outputimages[ind])
  # Run the proposed method (PEM) on SLURM (Multi-atlas segmentation)
  #sbatch(cmd_pem, mem=slurm_memory, n=slurm_nthreads, c=slurm_ncores, queue=slurm_queue, verbose=True, dryrun=False)

  if not (os.path.exists(outputimages[ind])):
      # run it on computer labs
      cmd_pem   = '/usr/lib/matlab/R2014a/bin/matlab -nodesktop -nosplash -r \\"cd(\'\\\'\'{0}\'\\\'\'); edgesTestingDemo(\'\\\'\'{1}\'\\\'\',\'\\\'\'{2}\'\\\'\',\'\\\'\'{3}\'\\\'\'); quit;\\"'.format(source_dir,modelname,inputimages[ind],output_dir)
      os.system('pmk.py --start --command "{0}"'.format(cmd_pem))