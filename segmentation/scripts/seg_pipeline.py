#!/usr/bin/python

"""
Cardiac segmentation pipeline

Author: Wenjia Bai
First created: 2015.03.31
Last modified: 2015.07.02 by wbai
"""

import os
import cmaseg

# Enter the data directory
data_root = '/vol/vipdata/data/biobank/cardiac/nifti/'
os.chdir(data_root)
list = sorted(os.listdir(data_root))
#list = list[4:]

for data in list:
    print data

    # The atlas directory
    pipeline = cmaseg.pipeline()
    pipeline.irtk_dir = '/vol/medic02/users/wbai/git/irtk_andreas_bin/bin'
    pipeline.data_dir = os.path.join(data_root, data)
    pipeline.atlas_dir = '/vol/medic02/users/wbai/data/3d_multi_atlas/old_version'
    pipeline.template_dir = '/vol/medic02/users/wbai/data/3d_atlas/v4/nreg_intensity_iter/ref3'
    pipeline.n_atlas = 20
    pipeline.n_proc = 20

    pipeline.segmentation()
