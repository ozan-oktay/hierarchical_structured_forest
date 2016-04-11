#!/usr/bin/python

import os
import tempfile as tf

def mkdirfun(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
def callmyfunction(mycmd):
    import subprocess
    cmd = subprocess.Popen(mycmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE)
    stdoutput = cmd.communicate()[0].strip('\n')
    print stdoutput
    return stdoutput
def removeContent(directory):
    import glob
    dirfiles = glob.glob(directory+'/*')
    for f in dirfiles:
        if os.path.isfile(f):
            os.remove(f)
def list2string(input_list):
    output = ''
    for element in input_list:
        output += ' {0}'.format(element)
    return output
def convertVol2Seq(inputnifti, outputnifti):
    import tempfile
    import shutil
    tmpdirpath_frames = tempfile.mkdtemp()
    callmyfunction('splitvolume {0} {1}/frame -slice'.format(inputnifti, tmpdirpath_frames))
    callmyfunction('makesequence {0}/frame* {1}'.format(tmpdirpath_frames, outputnifti))
    shutil.rmtree(tmpdirpath_frames)
def numberOfFiles(directory):
    path, dirs, files = os.walk(directory).next()
    return len(files)
def copyFile(src, dest):
    import shutil
    try:
        shutil.copy(src, dest)
    # eg. src and dest are the same file
    except shutil.Error as e:
        print('Error: %s' % e)
    # eg. source or destination doesn't exist
    except IOError as e:
        print('Error: %s' % e.strerror)

# Input directories
source_directory = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mritrainingdata'
target_directory = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mritrainingdata_pre_aligned'; mkdirfun(target_directory)

# Get the list of the images
source_img_dir   = os.path.join(source_directory,'images')
source_img_names = sorted(os.walk(source_img_dir).next()[2])

# Reference Image / Copy it to the target space
source_ref_dir = os.path.join(source_directory,'reference')
target_ref_dir = os.path.join(target_directory,'reference'); mkdirfun(target_ref_dir)
callmyfunction('cp -rf {0} {1}'.format(source_ref_dir,target_directory))
reference_image = os.path.join(target_ref_dir,'ref_image.nii.gz')

# Target Directories
target_img_dir = os.path.join(target_directory,'images');      mkdirfun(target_img_dir)
target_bou_dir = os.path.join(target_directory,'boundary');    mkdirfun(target_bou_dir)
target_grn_dir = os.path.join(target_directory,'groundtruth'); mkdirfun(target_grn_dir)
target_lm_dir  = os.path.join(target_directory,'landmarks');   mkdirfun(target_lm_dir)
target_dof_dir = os.path.join(target_directory,'dof_to_ref');  mkdirfun(target_dof_dir)

# Loop over all the source images
for source_img_name in source_img_names:

    # Print source image
    patient_name = source_img_name.split('.nii.gz')[0]; print source_img_name
    source_img_name = os.path.join(source_img_dir,source_img_name)

    # Perform Sampling to the same image resolution
    tf1 = tf.NamedTemporaryFile(mode = 'w+b', suffix = '.nii.gz', prefix = 'tmp', delete = True); tmp_target_img_name = tf1.name
    callmyfunction('transformation {0} {1} -linear -target {2} -matchInputType'.format(source_img_name,tmp_target_img_name,reference_image))

    # Robust Block Matching Based registration
    tf2 = tf.NamedTemporaryFile(mode = 'w+b', suffix = '.txt', prefix = 'tmp', delete = True); tmp_rbm_txt_name = tf2.name
    target_dof_name = os.path.join(target_dof_dir,patient_name+'.dof.gz')
    callmyfunction('reg_aladin -ref {0} -flo {1} -ln 3 -maxit 15 -nac -aff {2} -omp 64'.format(reference_image,tmp_target_img_name,tmp_rbm_txt_name))

    # Generate the dof
    callmyfunction('niftk2dof {0} {1}'.format(tmp_rbm_txt_name,target_dof_name))

    # Transform the boundary map - groundtruth - images - landmarks
    target_img_name = os.path.join(target_img_dir,patient_name+'.nii.gz')
    target_bou_name = os.path.join(target_bou_dir,patient_name+'.nii.gz')
    target_grn_name = os.path.join(target_grn_dir,patient_name+'.nii.gz')
    target_lm_name  = os.path.join(target_lm_dir,patient_name+'.vtk')

    source_bou_name = os.path.join(os.path.join(source_directory,'boundary'),patient_name+'.nii.gz')
    source_grn_name = os.path.join(os.path.join(source_directory,'groundtruth'),patient_name+'.nii.gz')
    source_lm_name  = os.path.join(os.path.join(source_directory,'landmarks'),patient_name+'.vtk')

    callmyfunction('transformation {0} {1} -linear -dofin {2} -matchInputType'.format(source_img_name,target_img_name,target_dof_name))
    callmyfunction('transformation {0} {1} -linear -dofin {2} -matchInputType'.format(source_bou_name,target_bou_name,target_dof_name))
    callmyfunction('transformation {0} {1} -linear -dofin {2} -matchInputType'.format(source_grn_name,target_grn_name,target_dof_name))
    callmyfunction('ptransformation {0} {1} -dofin {2} -invert'.format(source_lm_name,target_lm_name,target_dof_name))

    # Create a temporary folder for the inference
    target_tmp_dir = os.path.join(target_directory,'tmp'); mkdirfun(target_tmp_dir)