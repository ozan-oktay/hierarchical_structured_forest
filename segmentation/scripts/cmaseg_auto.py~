# Cardiac segmentation pipeline
import os, multiprocessing, re

#
# Indirect call for multi-processing of instance methods
#
def call_it(instance, name, args):
    # If we directly call the instance method inside the class, the following error would occur:
    # "PicklingError: Can't pickle <type 'instancemethod'>: attribute lookup __builtin__.instancemethod failed"
    #
    # Disadvantage: each time, an instance will be copied and destructed. So the destructor will be called multiple times.
    # http://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-pythons-multiprocessing-pool-ma
    return getattr(instance, name)(*args)


#
# Class representing the segmentation pipeline
#
class pipeline(object):
    
    def __init__(self, parent = None):
        super(pipeline, self).__init__()

    def reg(self, i, prev_method, cur_method, fr, mask): 
        if re.match('p\w*', cur_method):
            # Point-based registration
            target = '{0}/landmarks.vtk'.format(self.data_dir)
            atlas = '{0}/Atlas{1:02d}/landmarks.vtk'.format(self.atlas_dir, i)
        elif re.match('\w*_label', cur_method):
            # Label-based registration
            target = '{0}/fusion_pbcc_nreg_{1}.nii.gz'.format(self.results_dir, fr)
            atlas = '{0}/Atlas{1:02d}/segmentation_{2}.nii.gz'.format(self.atlas_dir, i, fr)
        else:
            # Image-based non-rigid registration
            target = '{0}/cine_sax_{1}_crop.nii.gz'.format(self.results_dir, fr)
            atlas = '{0}/Atlas{1:02d}/lvsa_{2}.nii.gz'.format(self.atlas_dir, i, fr)

        method = re.sub('_label', '', cur_method)
        cur_dof = '{0}/{1}_to_Atlas{2:02d}_{3}.dof.gz'.format(self.dof_dir, cur_method, i, fr)

        # prreg does only have a ED dof file
        if re.match('p\w', prev_method):
            prev_dof = '{0}/{1}_to_Atlas{2:02d}_ED.dof.gz'.format(self.dof_dir, prev_method, i)
        else:
            prev_dof = '{0}/{1}_to_Atlas{2:02d}_{3}.dof.gz'.format(self.dof_dir, prev_method, i, fr)

        # prreg does not require the parameter file
        if re.match('p\w', cur_method):
            par = ''
        else:
            par = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'par/{0}.cfg'.format(cur_method))

        # Use mask for acceleration
        if mask:
            mask = '{0}/{1}.nii.gz'.format(self.results_dir, mask)

        par_str = ''
        if par:
            par_str += '-parin {0} '.format(par)
        if prev_method:
            par_str += '-dofin {0} '.format(prev_dof)
        if mask:
            par_str += '-mask {0} '.format(mask)

        os.system('{0}/{1} {2} {3} -dofout {4} {5}'.format(self.irtk_dir, method, target, atlas, cur_dof, par_str))

    def transform(self, i, method, fr, crop):
        # Image transformation
        if crop:
            target_image = '{0}/cine_sax_{1}_{2}.nii.gz'.format(self.results_dir, fr, crop)
        else:
            target_image = '{0}/cine_sax_{1}.nii.gz'.format(self.results_dir, fr)
        atlas_image = '{0}/Atlas{1:02d}/lvsa_{2}.nii.gz'.format(self.atlas_dir, i, fr)
        atlas_label = '{0}/Atlas{1:02d}/segmentation_{2}.nii.gz'.format(self.atlas_dir, i, fr)
        
        dof = '{0}/{1}_to_Atlas{2:02d}_{3}.dof.gz'.format(self.dof_dir, method, i, fr)
        image_prop = '{0}/image_from_Atlas{1:02d}_{2}_{3}.nii.gz'.format(self.warped_dir, i, method, fr)
        os.system('{0}/transformation {1} {2} -linear -dofin {3} -target {4}'.format(self.irtk_dir, atlas_image, image_prop, dof, target_image))
        label_prop = '{0}/label_from_Atlas{1:02d}_{2}_{3}.nii.gz'.format(self.warped_dir, i, method, fr)
        os.system('{0}/transformation {1} {2} -nn -dofin {3} -target {4}'.format(self.irtk_dir, atlas_label, label_prop, dof, target_image))

    def segmentation(self):
        # Create result directory
        self.results_dir = '{0}/results'.format(self.data_dir)
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)

        self.dof_dir = '{0}/results/dof'.format(self.data_dir)
        if not os.path.exists(self.dof_dir):
            os.mkdir(self.dof_dir)

        self.warped_dir = '{0}/results/warped'.format(self.data_dir)
        if not os.path.exists(self.warped_dir):
            os.mkdir(self.warped_dir)

        # The first frame is ED frame
        ED = 0

        # Detect ES frame
        command_return = os.popen('detect_ES_frame {0}/cine_sax.nii.gz'.format(self.data_dir))
        ES = int(command_return.read())
        command_return.close()

        # Split sequence
        os.system('{0}/splitvolume {1}/cine_sax.nii.gz {2}/cine_sax_ -sequence'.format(self.irtk_dir, self.data_dir, self.results_dir))

        # Process ED and ES frames
        os.system('cp {0}/cine_sax_{1:02d}.nii.gz {0}/cine_sax_ED.nii.gz'.format(self.results_dir, ED))
        os.system('cp {0}/cine_sax_{1:02d}.nii.gz {0}/cine_sax_ES.nii.gz'.format(self.results_dir, ES))
        os.system('rm -f {0}/cine_sax_[0-9]*.nii.gz'.format(self.results_dir))

        # Select landmarks
        #print 'Please select 6 landmarks at the ED phase'
        #os.system('{0}/rview {1}/cine_sax_00.nii.gz'.format(self.irtk_dir, self.results_dir))

        # Check landmarks

        # Pool for parallel processing
        pool = multiprocessing.Pool(self.n_proc)

        #
        # This step may be replaced by image-based registration using the original orientation information
        #
        # Landmark registration
        results = [pool.apply_async(call_it, (self, 'reg', (i,'','prreg','ED',''))) for i in range(1, self.n_atlas+1)]
        for r in results:
            r.get()
        
        # Image transformation based on landmark registration
        results = [pool.apply_async(call_it, (self, 'transform', (i,'prreg','ED',''))) for i in range(1, self.n_atlas+1)]
        for r in results:
            r.get()
        
        # Initial segmentation to detect the ROI
        # Transformed atlas images and label maps
        image_props = ''
        label_props = ''
        for i in range(1, self.n_atlas+1):
            image_prop = '{0}/image_from_Atlas{1:02d}_prreg_ED.nii.gz'.format(self.warped_dir, i)
            image_props += image_prop + ' '
            label_prop = '{0}/label_from_Atlas{1:02d}_prreg_ED.nii.gz'.format(self.warped_dir, i)
            label_props += label_prop + ' '
        
        # Label fusion
        target_image = '{0}/cine_sax_ED.nii.gz'.format(self.results_dir)
        fusion = '{0}/fusion_mv_prreg_ED.nii.gz'.format(self.results_dir)
        os.system('label_fusion {0} {1} {2} {3} {4} -method MV'.format(target_image, self.n_atlas, image_props, label_props, fusion))
        
        # Crop the images to save computation for image registration
        fusion_crop = '{0}/fusion_mv_prreg_ED_crop.nii.gz'.format(self.results_dir)
        os.system('auto_crop_image {0} {1} -reserve 20'.format(fusion, fusion_crop))

        target_image = '{0}/cine_sax_ED.nii.gz'.format(self.results_dir)
        target_image_crop = '{0}/cine_sax_ED_crop.nii.gz'.format(self.results_dir)
        os.system('{0}/region {1} {2} -ref {3}'.format(self.irtk_dir, target_image, target_image_crop, fusion_crop))

        target_image = '{0}/cine_sax_ES.nii.gz'.format(self.results_dir)
        target_image_crop = '{0}/cine_sax_ES_crop.nii.gz'.format(self.results_dir)
        os.system('{0}/region {1} {2} -ref {3}'.format(self.irtk_dir, target_image, target_image_crop, fusion_crop))
        
        # Mask for the heart region to save computation and improve accuracy
        mask = '{0}/mask.nii.gz'.format(self.results_dir)
        os.system('{0}/padding {1} {1} {2} 0 1 -invert'.format(self.irtk_dir, fusion_crop, mask))
        os.system('{0}/dilation {1} {1} -iterations 7'.format(self.irtk_dir, mask, mask))
        
        # Segmentation at ED and ES
        for fr in ['ED', 'ES']:
            # Affine registration
            results = [pool.apply_async(call_it, (self, 'reg', (i,'prreg','areg',fr,''))) for i in range(1, self.n_atlas+1)]
            for r in results:
                r.get()
            
            # Non-rigid registration
            results = [pool.apply_async(call_it, (self, 'reg', (i,'areg','nreg',fr,'mask'))) for i in range(1, self.n_atlas+1)]
            for r in results:
               r.get()
            
            # Image transformation
            results = [pool.apply_async(call_it, (self, 'transform', (i,'nreg',fr,'crop'))) for i in range(1, self.n_atlas+1)]
            for r in results:
               r.get()
            
            # TODO: replace nreg by ireg
            
            # Transformed atlas images and label maps
            image_props = ''
            label_props = ''
            for i in range(1, self.n_atlas+1):
                image_prop = '{0}/image_from_Atlas{1:02d}_nreg_{2}.nii.gz'.format(self.warped_dir, i, fr)
                image_props += image_prop + ' '
                label_prop = '{0}/label_from_Atlas{1:02d}_nreg_{2}.nii.gz'.format(self.warped_dir, i, fr)
                label_props += label_prop + ' '
            
            # Label fusion
            target_image = '{0}/cine_sax_{1}_crop.nii.gz'.format(self.results_dir, fr)
            fusion = '{0}/fusion_pbcc_nreg_{1}.nii.gz'.format(self.results_dir, fr)
            par = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'par/pbcc.cfg')
            os.system('label_fusion {0} {1} {2} {3} {4} -method PBCC -par {5}'.format(target_image, self.n_atlas, image_props, label_props, fusion, par))
            
            # Refine registration using label maps
            # Affine registration
            results = [pool.apply_async(call_it, (self, 'reg', (i,'prreg','areg_label',fr,''))) for i in range(1, self.n_atlas+1)]
            for r in results:
                r.get()
            
            # Non-rigid registration
            results = [pool.apply_async(call_it, (self, 'reg', (i,'areg_label','nreg_label',fr,''))) for i in range(1, self.n_atlas+1)]
            for r in results:
                r.get()
            
            # Image transformation
            results = [pool.apply_async(call_it, (self, 'transform', (i,'nreg_label',fr,'crop'))) for i in range(1, self.n_atlas+1)]
            for r in results:
                r.get()
            
            # Transformed atlas images and label maps
            image_props = ''
            label_props = ''
            for i in range(1, self.n_atlas+1):
                image_prop = '{0}/image_from_Atlas{1:02d}_nreg_label_{2}.nii.gz'.format(self.warped_dir, i, fr)
                image_props += image_prop + ' '
                label_prop = '{0}/label_from_Atlas{1:02d}_nreg_label_{2}.nii.gz'.format(self.warped_dir, i, fr)
                label_props += label_prop + ' '
            
            # Label fusion
            target_image = '{0}/cine_sax_{1}_crop.nii.gz'.format(self.results_dir, fr)
            fusion = '{0}/fusion_pbcc_nreg_label_{1}.nii.gz'.format(self.results_dir, fr)
            par = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'par/pbcc.cfg')
            os.system('label_fusion {0} {1} {2} {3} {4} -method PBCC -par {5}'.format(target_image, self.n_atlas, image_props, label_props, fusion, par))

            # Fit to the template
            target = '{0}/landmarks.vtk'.format(self.data_dir)
            template = '{0}/landmarks_in_ref3_nreg.vtk'.format(self.template_dir)
            prreg_dof = '{0}/prreg_to_template.dof.gz'.format(self.dof_dir)
            os.system('{0}/prreg {1} {2} -dofout {3}'.format(self.irtk_dir, target, template, prreg_dof))

            target_label = fusion
            template_label = '{0}/segmentation_{1}_in_ref3_nreg.nii.gz'.format(self.template_dir, fr)
            areg_dof = '{0}/areg_label_to_template_{1}.dof.gz'.format(self.dof_dir, fr)
            par = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'par/areg_label.cfg')
            mask = '{0}/mask.nii.gz'.format(self.results_dir)
            os.system('{0}/areg {1} {2} -parin {3} -dofin {4} -dofout {5}'.format(self.irtk_dir, target_label, template_label, par, prreg_dof, areg_dof))

            # Refine areg by non-rigid surface registration which can warp the mesh model more accurately
            myo = '{0}/myo_pbcc_nreg_{1}.nii.gz'.format(self.results_dir, fr)
            os.system('{0}/padding {1} {1} {2} 2 0 -invert'.format(self.irtk_dir, fusion, myo))        

            myo_mesh = '{0}/myo_pbcc_nreg_{1}.vtk'.format(self.results_dir, fr)
            os.system('{0}/mcubes {1} {2} 1'.format(self.irtk_dir, myo, myo_mesh))

            template_myo_mesh = '{0}/myo_{1}_sm.vtk'.format(self.template_dir, fr)
            snreg_dof = '{0}/snreg_to_template_{1}.dof.gz'.format(self.dof_dir, fr)
            os.system('{0}/snreg {1} {2} -dofin {3} -dofout {4} -locator 0 -ds 5'.format(self.irtk_dir, myo_mesh, template_myo_mesh, areg_dof, snreg_dof))

            # Transformation of template label map, mesh and AHA model
            template_label = '{0}/segmentation_{1}_in_ref3_nreg.nii.gz'.format(self.template_dir, fr)
            label_prop = '{0}/template_fit_snreg_{1}.nii.gz'.format(self.results_dir, fr)
            os.system('{0}/transformation {1} {2} -nn -dofin {3} -target {4}'.format(self.irtk_dir, template_label, label_prop, snreg_dof, target_label))

            template_mesh = '{0}/heart_{1}_sm.vtk'.format(self.template_dir, fr)
            mesh_prop = '{0}/template_fit_snreg_{1}.vtk'.format(self.results_dir, fr)
            os.system('{0}/ptransformation {1} {2} -dofin {3} -invert'.format(self.irtk_dir, template_mesh, mesh_prop, snreg_dof))
                
            if fr == 'ED':
                # AHA17 model is only available at ED
                aha17_model = '{0}/AHA17_segment_model.nii.gz'.format(self.template_dir)
                model_prop = '{0}/AHA17_fit_snreg_{1}.nii.gz'.format(self.results_dir, fr)
                os.system('{0}/transformation {1} {2} -nn -dofin {3} -target {4}'.format(self.irtk_dir, aha17_model, model_prop, snreg_dof, target_label))
            
                aha17_mesh = '{0}/myo_{1}_sm_AHA17.vtk'.format(self.template_dir, fr)
                mesh_prop = '{0}/AHA17_fit_snreg_{1}.vtk'.format(self.results_dir, fr)
                os.system('{0}/ptransformation {1} {2} -dofin {3} -invert'.format(self.irtk_dir, aha17_mesh, mesh_prop, snreg_dof))
            
        # Close the pool
        pool.close()
        pool.join()
        
