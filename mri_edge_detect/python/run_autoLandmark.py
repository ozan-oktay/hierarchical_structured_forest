__author__ = 'oo2113'
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

def get_filepaths(directory,strendswith=''):
    """
    This function will generate the file names in a directory
    tree by walking the tree either top-down or bottom-up. For each
    directory in the tree rooted at directory top (including top itself),
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            if filename.endswith(strendswith):
                # Join the two strings in order to form the full filepath.
                filepath = os.path.join(root, filename)
                file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.

# Algorithm Parameters
firstClassifier = 'mriFirst_hier'
secondClassifier= 'mriSecond_hier_E4'
slurm_ncores    = 8
slurm_nthreads  = 1
slurm_memory    = 25
slurm_queue     = 'short'

# Directories
target_dir   = '/vol/biomedic/users/oo2113/str_hier_forest_mri/segmentation/targets'
source_dir   = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mri_edge_detect'
inputimages  = get_filepaths(target_dir,'cine_sax_ED.nii.gz')

for inputimage in inputimages:

    pemOut1   = os.path.dirname(inputimage) + '/pems';  createFolder(pemOut1)
    pemOut2   = os.path.dirname(inputimage) + '/pems2'; createFolder(pemOut2)

    if not os.path.exists(pemOut2+'/cine_sax_ED_pem.nii.gz'):
        cmd_pem1  = '/usr/lib/matlab/R2014a/bin/matlab -nodesktop -nosplash -r \\"cd(\'{0}\'); addpath(\'{0}\'); edgesTestingDemo(\'{1}\',\'{2}\',\'{3}\'); quit;\\"'.format(source_dir,firstClassifier, inputimage,pemOut1)
        cmd_pem2  = '/usr/lib/matlab/R2014a/bin/matlab -nodesktop -nosplash -r \\"cd(\'{0}\'); addpath(\'{0}\'); edgesTestingDemo(\'{1}\',\'{2}\',\'{3}\'); quit;\\"'.format(source_dir,secondClassifier,inputimage,pemOut2)
        cmd_slurm = cmd_pem1+'; '+cmd_pem2
        sbatch(cmd_slurm, mem=slurm_memory, n=slurm_nthreads, c=slurm_ncores, queue=slurm_queue, verbose=True, dryrun=False)
