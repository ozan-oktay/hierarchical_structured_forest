#!/usr/bin/python

"""
Report cardiac function in a spreadsheet

Author: Wenjia Bai
First created: 2015.05.06
Last modified: 2015.05.06 by wbai
"""

import os, csv, string

# Enter the data directory
data_root = '/vol/vipdata/data/biobank/cardiac/nifti/'
os.chdir(data_root)

# Read the names
f = open('../names.csv', 'r')
reader = csv.reader(f)

# Write the EDV, ESV, LVM, etc
f2 = open('../cardiac_report.csv', 'w')
writer = csv.writer(f2)
row2 = ['Patient ID', 'Name', 'LVEDV (mL)', 'LVESV (mL)', 'LVEJ (%)', 'LVM (g)', 'RVEDV (mL)', 'RVESV (mL)', 'RVEJ (%)']
writer.writerow(row2)

for row in reader:
    if row == ['Patient ID', 'Name']:
        continue
    
    data = row[1]
    data = string.replace(data, ' ', '') + '_'
    print data

    data_dir = '{0}/results'.format(data)
    
    # ED
    seg = '{0}/fusion_pbcc_nreg_label_ED.nii.gz'.format(data_dir)
    command_return = os.popen('measure_volume {0} 1'.format(seg))
    lv_edv = float(command_return.read())
    command_return.close()

    command_return = os.popen('measure_volume {0} 2'.format(seg))
    myo_density = 1.05
    lvm = float(command_return.read()) * myo_density
    command_return.close()

    command_return = os.popen('measure_volume {0} 4'.format(seg))
    rv_edv = float(command_return.read())
    command_return.close()

    # ES
    seg = '{0}/fusion_pbcc_nreg_label_ES.nii.gz'.format(data_dir)
    command_return = os.popen('measure_volume {0} 1'.format(seg))
    lv_esv = float(command_return.read())
    command_return.close()

    command_return = os.popen('measure_volume {0} 4'.format(seg))
    rv_esv = float(command_return.read())
    command_return.close()

    # Ejection fraction
    lv_ej = '{0:.3f}'.format((lv_edv - lv_esv) / lv_edv * 100)
    rv_ej = '{0:.3f}'.format((rv_edv - rv_esv) / rv_edv * 100)

    # Write to the csv file
    row2 = row + [lv_edv, lv_esv, lv_ej, lvm, rv_edv, rv_esv, rv_ej]
    writer.writerow(row2)
