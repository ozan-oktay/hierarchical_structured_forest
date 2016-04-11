__author__ = 'oo2113'
import os
import numpy as np

def callmyfunction(mycmd):
    import subprocess
    cmd = subprocess.Popen(mycmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE)
    stdoutput = cmd.communicate()[0].strip('\n')
    print stdoutput
    return stdoutput
def mkdirfun(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

# Input directories
dataset_dir      = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mritestingdata'
segmentation_dir = dataset_dir + '/groundtruth'
landmark_dir     = dataset_dir + '/landmarks'
reference_lm     = dataset_dir + '/reference/ref_lm.vtk'
reference_img    = dataset_dir + '/reference/ref_image_cropped.nii.gz'
tmp_dir          = '/vol/bitbucket/oo2113/tmp_lle'; mkdirfun(tmp_dir)

# Retrieve the list of segmentation images and landmark points
segmentation_names = sorted(os.walk(segmentation_dir).next()[2]); segmentation_names = [name for name in segmentation_names if name.endswith('.nii.gz')]
landmark_names     = sorted(os.walk(landmark_dir).next()[2]);     landmark_names     = [name for name in landmark_names     if name.endswith('.vtk')]

###### Remove this part of the code ######
segmentation_names = segmentation_names[:10]; print segmentation_names
landmark_names     = landmark_names[:10];     print landmark_names
###### ###### ###### ###### ###### ######

assert ( len(segmentation_names) == len(landmark_names) )
num_subjects       = len(segmentation_names)

# Initialize the proximity matrix for all the images
proximity_matrix = np.zeros((num_subjects,num_subjects))

# Translate the images and store them in a tmp dir
for segmentation_name, lm_name in zip(segmentation_names,landmark_names):
    assert (segmentation_name.split('_')[0] == lm_name.split('_')[0])
    tmp_segmentation_name = tmp_dir + '/' + segmentation_name
    tmp_dof_name          = tmp_dir + '/' + segmentation_name.split('.nii.gz')[0] + '.dof.gz'
    callmyfunction('prreg {0} {1} -dofout {2}'.format(reference_lm,landmark_dir+'/'+lm_name,tmp_dof_name))
    callmyfunction('transformation {0} {1} -dofin {2} -nn -matchInputType -target {3}'.format(segmentation_dir+'/'+segmentation_name,tmp_segmentation_name,tmp_dof_name,reference_img))

# Compute the dice metric / Hausdorff distance between the segmentation pairs
tmp_segmentation_names = sorted(os.walk(tmp_dir).next()[2]); tmp_segmentation_names = [tmp_dir+'/'+name for name in tmp_segmentation_names if name.endswith('.nii.gz')]
assert ( len(tmp_segmentation_names) == len(segmentation_names) )
for first_ind,target_img in enumerate(tmp_segmentation_names):
    for second_ind,source_img in enumerate(tmp_segmentation_names):

        # Dice computation
        avg_dice = 0.0; label_ids = [0,1,4]
        for label_id in label_ids:
            cmd_out    = callmyfunction('dicemetric {0} {1} -minvalue {2} -maxvalue {2}'.format(target_img,source_img,label_id))
            dice_score = (float(cmd_out.split('Dice metric of frame 0 is ')[1]))
            avg_dice  += ( dice_score / float(len(label_ids)) )

        # Populate the proximity matrix and perform laplacian eigenmaps
        proximity_matrix[first_ind,second_ind] = avg_dice

print proximity_matrix

# Perform manifold learning
import sklearn.manifold.SpectralEmbedding

#spec = SpectralEmbedding(n_components=2, affinity='precomputed', gamma=None, random_state=None, eigen_solver=None, n_neighbors=None)[source]
#spec.fit ()
#fit(X, y=None)

