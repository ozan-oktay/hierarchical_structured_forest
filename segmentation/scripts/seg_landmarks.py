#!/usr/bin/python

"""
Cardiac landmarks selection (manual process)

Author: Wenjia Bai
First created: 2015.03.31
Last modified: 2015.04.02 by wbai
"""
import os

# Enter the data directory
data_root = '/vol/vipdata/data/biobank/cardiac/nifti/'
os.chdir(data_root)

list = sorted(os.listdir(data_root))
list = list[2:]

for data in list:
    print data
    os.chdir(data)
    os.system('rview cine_sax.nii.gz')
    os.system('mv l.vtk landmarks.vtk')
    os.chdir('..')
