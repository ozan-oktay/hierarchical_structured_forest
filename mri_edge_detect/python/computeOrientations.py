#!/usr/bin/python

__author__ = 'ozan'
import os
import subprocess
import irtk
import scipy.io
import argparse

# Function call
def callmyfunction(mycmd):
    subprocess.call(mycmd, shell=True, executable="/bin/bash")

# Create directory function
def mycreatedir (l_directory):
    if not os.path.exists(l_directory):
        os.makedirs(l_directory)

parser   = argparse.ArgumentParser(description='')
parser.add_argument( '--indir',  type=str, help='',  required=False )
args     = parser.parse_args()
inputdir = args.indir

if not inputdir:
    inputdir = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mritestingdata'

mov_landmarks  = []
dofoutnames    = []
z_rotations    = []
filenames      = []
ref_landmark     = inputdir + '/reference/ref_lm.vtk'
mov_landmark_dir = inputdir + '/landmarks'
out_dof_dir      = inputdir + '/dofs'; mycreatedir(out_dof_dir)

for dirpath, subdirs, files in os.walk(mov_landmark_dir):
    for l_file in files:
        if l_file.endswith('vtk'):
            mov_landmarks.append(dirpath+'/'+l_file)
            dofoutnames.append(out_dof_dir+'/'+l_file.split('.vtk')[0]+'.dof.gz')
            filenames.append(dirpath+'/'+l_file.split('.vtk')[0]+'_lm.vtk')

mov_landmarks = sorted(mov_landmarks)
dofoutnames   = sorted(dofoutnames)
filenames     = sorted(filenames)

for mov_landmark,dofoutname in zip(mov_landmarks,dofoutnames):
    callmyfunction('prreg {0} {1} -dofout {2}'.format(ref_landmark,mov_landmark,dofoutname))
    dofparams=irtk.AffineTransformation(filename=dofoutname)
    print dofparams
    z_rotations.append(dofparams.rz)

matdict={'filename':filenames, 'z_rot':z_rotations}
scipy.io.savemat(out_dof_dir+'/rotation.mat', matdict )