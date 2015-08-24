#!/usr/bin/python

"""
Cardiac motion tracking pipeline

Author: Wenjia Bai
First created: 2015.03.31
Last modified: 2015.04.02 by wbai
"""

import os

# Enter the data directory
data_root = '/vol/vipdata/data/biobank/cardiac/nifti/'
os.chdir(data_root)
list = sorted(os.listdir(data_root))
list = list[10:]

irtk_dir = '/vol/medic02/users/wbai/git/irtk_andreas_bin/bin'
par = '~/Work/Programming/python/biobank/par/motion_ffd4d.cfg'

for data in list:
    print data
    data_dir = os.path.join(data_root, data)

    # Create result directory
    results_dir = '{0}/results'.format(data_dir)
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    track_dir = '{0}/results/track'.format(data_dir)
    if not os.path.exists(track_dir):
        os.mkdir(track_dir)

    # Crop the image sequence
    seq = '{0}/cine_sax.nii.gz'.format(data_dir)
    seq_crop = '{0}/cine_sax_crop.nii.gz'.format(track_dir)
    image_crop = '{0}/cine_sax_ED_crop.nii.gz'.format(results_dir)
    os.system('{0}/region {1} {2} -ref {3}'.format(irtk_dir, seq, seq_crop, image_crop))

    # Split sequence
    os.system('{0}/splitvolume {1}/cine_sax_crop.nii.gz {1}/cine_sax_crop_ -sequence'.format(irtk_dir, track_dir))

    # Count number of frames
    command_return = os.popen('ls {0}/cine_sax_crop_[0-9]*.nii.gz | wc -l'.format(track_dir))
    num_fr = int(command_return.read())
    command_return.close()

    # Perform motion tracking
    images  = ''
    for fr in range(0, num_fr):
        image   = '{0}/cine_sax_crop_{1:02d}.nii.gz'.format(track_dir, fr)
        images += image + ' '
    images += '{0}/cine_sax_crop_00.nii.gz'.format(track_dir)

    dof = '{0}/motion_ffd4d.dof.gz'.format(results_dir)
    log = '{0}/motion_ffd4d.log'.format(results_dir)
    os.system('motion_track {0} {1} {2} -parin {3} -gpu on -display_time -display_opt > {4}'.format(num_fr + 1, images, dof, par, log))
    if os.path.exists('{0}.gz'.format(log)):
        os.system('rm -f {0}.gz'.format(log))
    os.system('gzip {0}'.format(log))

    # Propagate surface mesh
    mesh_ED = '{0}/template_fit_snreg_ED.vtk'.format(results_dir)

    for fr in range(0, num_fr):
        mesh_out = '{0}/template_fit_snreg_prop_{1:02d}.vtk'.format(track_dir, fr)
        if fr == 0:
            os.system('cp {0} {1}'.format(mesh_ED, mesh_out))
        else:
            os.system('ptransformation {0} {1} -dofin {2} -time {3}'.format(mesh_ED, mesh_out, dof, fr))
