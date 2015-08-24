# Cardiac segmentation report
import os, time
from reportlab.pdfgen import canvas

def report(data_name, data_dir, output):
    print data_dir
    print output

    # ED
    seg = '{0}/fusion_pbcc_nreg_label_ED.nii.gz'.format(data_dir)
    command_return = os.popen('measure_volume {0} 1'.format(seg))
    lv_edv = float(command_return.read())
    command_return.close()

    command_return = os.popen('measure_volume {0} 2'.format(seg))
    myo_density = 1.05
    lvm_ed = float(command_return.read()) * myo_density
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
    lv_ej = (lv_edv - lv_esv) / lv_edv * 100
    rv_ej = (rv_edv - rv_esv) / rv_edv * 100

    # Segmentation screenshot
    irtk_dir = '/vol/medic02/users/wbai/git/irtk_andreas_bin/bin'
    os.system('{0}/display {1}/cine_sax_ED.nii.gz -seg {1}/fusion_pbcc_nreg_label_ED.nii.gz -lut /homes/wbai/Work/Programming/python/biobank/par/lut.cfg -offscreen {1}/seg_ED.png'.format(irtk_dir, data_dir))
    os.system('{0}/display {1}/cine_sax_ES.nii.gz -seg {1}/fusion_pbcc_nreg_label_ES.nii.gz -lut /homes/wbai/Work/Programming/python/biobank/par/lut.cfg -offscreen {1}/seg_ES.png'.format(irtk_dir, data_dir))

    os.system('{0}/display {1}/cine_sax_ED.nii.gz -object {1}/template_fit_snreg_ED.vtk -offscreen {1}/temp_fit_ED.png'.format(irtk_dir, data_dir))
    os.system('{0}/display {1}/cine_sax_ES.nii.gz -object {1}/template_fit_snreg_ES.vtk -offscreen {1}/temp_fit_ES.png'.format(irtk_dir, data_dir))

    # Drawing boundary
    left  = 50
    right = 550
    top   = 800

    # Start report
    c = canvas.Canvas(output)
    c.drawString(left, top+5, "Cardiac report")
    c.drawString(left+300, top+5, "Generated at {0}".format(time.strftime("%c")))
    c.line(left, top, right, top)
    c.drawString(left, top-15, 'Subject: {0}'.format(data_name))
     
    c.drawString(left+120, top-45, "End-diastole")
    c.drawString(left+200, top-45, "End-systole")
    c.drawString(left+280, top-45, "Ejec. frac.")

    c.drawString(left, top-60, " LV cavity (mL)")
    c.drawString(left+120, top-60, ' {0:.2f}'.format(lv_edv))
    c.drawString(left+200, top-60, ' {0:.2f}'.format(lv_esv))
    c.drawString(left+280, top-60, ' {0:.2f}%'.format(lv_ej))

    c.drawString(left, top-75, " LV myo (gram)")
    c.drawString(left+120, top-75, ' {0:.2f}'.format(lvm_ed))

    c.drawString(left, top-90, " RV cavity (mL)")
    c.drawString(left+120, top-90, ' {0:.2f}'.format(rv_edv))
    c.drawString(left+200, top-90, ' {0:.2f}'.format(rv_esv))
    c.drawString(left+280, top-90, ' {0:.2f}%'.format(rv_ej))
    
    img_width = (right - left) / 2.1
    c.drawInlineImage('{0}/seg_ED.png'.format(data_dir), left, top-120-img_width, img_width, img_width)
    c.drawInlineImage('{0}/seg_ES.png'.format(data_dir), left+img_width+10, top-120-img_width, img_width, img_width)
    c.drawString(left+50, top-135-img_width, '(a) Segmentation at ED')
    c.drawString(left+img_width+60, top-135-img_width, '(b) Segmentation at ES')

    c.drawInlineImage('{0}/temp_fit_ED.png'.format(data_dir), left, top-400-img_width, img_width, img_width)
    c.drawInlineImage('{0}/temp_fit_ES.png'.format(data_dir), left+img_width+10, top-400-img_width, img_width, img_width)
    c.drawString(left+50, top-415-img_width, '(a) Template fitting at ED')
    c.drawString(left+img_width+60, top-415-img_width, '(b) Template fitting at ES')

    c.showPage()
    c.save()
