#!/usr/bin/python

"""
Cardiac segmentation pipeline

Author: Wenjia Bai
First created: 2015.03.31
Last modified: 2015.08.03 by wbai

modified by: oo2113 - 2015.08.24
"""

import os
import cmaseg_auto

# Enter the data directory
data_root = '/vol/biomedic/users/oo2113/str_hier_forest_mri/segmentation/targets'
os.chdir(data_root)
list = sorted(os.listdir(data_root))

for loopId,data in enumerate(list):
    print data

    # The atlas directory
    pipeline = cmaseg_auto.pipeline()
    pipeline.irtk_dir = '/vol/medic02/users/wbai/git/irtk_andreas_bin/bin'
    pipeline.data_dir = os.path.join(data_root, data)
    pipeline.atlas_dir = '/vol/biomedic/users/oo2113/str_hier_forest_mri/segmentation/atlases'
    pipeline.template_dir = '/vol/medic02/users/wbai/data/3d_atlas/v4/nreg_intensity_iter/ref3'
    pipeline.targetId = (loopId+1) # leave-one out cross validation
    pipeline.n_atlas = 20
    pipeline.n_proc = 20

    pipeline.segmentation()
