#!/usr/bin/python

"""
Extract landmarks from Ozan's Hough vote map.
"""

import os

# Enter the data directory
data_root = '/vol/vipdata/data/biobank/cardiac/nifti/'
os.chdir(data_root)
list = sorted(os.listdir(data_root))

for data in list:
    print data

    hough_map = '/vol/biomedic/users/oo2113/for_Wenjia/{0}_ed0_3D_lm.nii.gz'.format(data)
    output = '{0}/landmarks_auto.vtk'.format(data)
    os.system('extract_landmarks {0} {1}'.format(hough_map, output))
    
