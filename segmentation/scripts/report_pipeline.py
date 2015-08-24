#!/usr/bin/python

"""
Cardiac report pipeline

Author: Wenjia Bai
First created: 2015.04.06
Last modified: 2015.06.02 by wbai
"""

import os
import cardiac_report

# Enter the data directory
data_root = '/vol/vipdata/data/biobank/cardiac/nifti/'
os.chdir(data_root)
list = sorted(os.listdir(data_root))
#list = list[:1]

pdfs = ''
for data in list:
    data_dir = '{0}/results'.format(data)
    output = '{0}/report.pdf'.format(data)
    cardiac_report.report(data, data_dir, output)
    pdfs += output + ' ' 
    
os.system('pdftk {0} cat output ../compiled_report.pdf'.format(pdfs))
