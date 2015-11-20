__author__ = 'oo2113'

import os
import irtk
from sklearn.ensemble import RandomForestRegressor
import numpy as np
import SimpleITK as sitk
import vtk
import scipy
import itertools

def vtk_point_read(filename):
    # load a vtk file as input
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    data   = reader.GetOutput()
    points = []
    for i in range(data.GetNumberOfPoints()):
        x, y, z = data.GetPoint(i)
        points.append([x,y,z])
    points = np.array(points,dtype=np.float64)
    return points

def comb2(n, k):
    from math import factorial
    return factorial(n) / factorial(k) / factorial(n - k)

# TRAINING PEM DIRECTORY
trainingImgDir   = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mritrainingdata_sec/images'
trainingPemDir   = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mritrainingdata_sec/pems'
trainingDofDir   = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mritrainingdata_sec/dofs'
trainingPemNames = []
trainingDofNames = []
trainingVtkNames = []
trainingImgNames = []

# TRAINING PARAMETERS
numAllSamples      = np.array([670])
numTrainingSamples = np.array([600])
inputSampleSize    = np.array([34,34,6])
downsampleSize     = 3
smoothWidth        = 2
numLandmarks       = 6

# INPUT TRAINING PEM IMAGES
for root,dirs,files in os.walk(trainingPemDir):
    for pemFile in files:
        if (pemFile.endswith('.nii.gz')) & ('lm' not in pemFile):
            trainingPemNames.append(root+'/'+pemFile)
        if (pemFile.endswith('.vtk')) & ('lm' in pemFile):
            trainingVtkNames.append(root+'/'+pemFile)
trainingPemNames = sorted(trainingPemNames)
trainingVtkNames = sorted(trainingVtkNames)

# INPUT TRAINING INTENSITY IMAGES
for root,dirs,files in os.walk(trainingImgDir):
    for imgFile in files:
        if imgFile.endswith('.nii.gz') & ('es' not in imgFile):
            trainingImgNames.append(root+'/'+imgFile)
trainingImgNames = sorted(trainingImgNames)

# ASSOCIATED DOF GROUND TRUTH INFORMATION
for root,dirs,files in os.walk(trainingDofDir):
    for dofFile in files:
        if dofFile.endswith('.dof.gz'):
            trainingDofNames.append(root+'/'+dofFile)
trainingDofNames = sorted(trainingDofNames)

# NUMBER OF TRAINING IMAGES
assert ( len(trainingPemNames) == len(trainingDofNames) )
assert ( len(trainingPemNames) == len(trainingVtkNames) )
assert ( len(trainingPemNames) == len(trainingImgNames) )
for trainingPemName,trainingDofName,trainingVtkName,trainingImgName in zip(trainingPemNames,trainingDofNames,trainingVtkNames,trainingImgNames):
    assert( trainingPemName.split('/')[-1].split('_pem.')[0] == trainingDofName.split('/')[-1].split('.')[0] )
    assert( trainingPemName.split('/')[-1].split('_pem.')[0] == trainingImgName.split('/')[-1].split('.')[0] )
    assert( trainingPemName.split('/')[-1].split('_pem.')[0] == trainingVtkName.split('/')[-1].split('.')[0].split('_lm')[0] )

numAllSamples  = np.min([numAllSamples,len(trainingPemNames)])
numPemFeatures = inputSampleSize[0]*inputSampleSize[1]*inputSampleSize[2]
numLMFeatures  = comb2(numLandmarks, 2) + comb2(comb2(numLandmarks, 2),2)
numHoGFeatures = 512 * inputSampleSize[2] * downsampleSize
numFeatures    = numLMFeatures + numPemFeatures + numHoGFeatures

collectedLabels = np.zeros([numAllSamples],dtype=float)
collectedFtrs   = np.zeros([numAllSamples,numFeatures],dtype=float)

# READ THE DOF FILES AND STORE GROUND TRUTH INFORMATION
for index in range(numAllSamples):

    dofname   = trainingDofNames[index]
    dofparams = irtk.AffineTransformation(filename=dofname)
    scale     = np.power(np.linalg.det(dofparams.matrix()),1.0/3.0)

    #collectedLabels[index] = scale
    collectedLabels[index] = dofparams.rz

    # READ THE GENERATED PEM FILE
    pemname = trainingPemNames[index]
    pemimg  = sitk.ReadImage(pemname)
    print 'processing image {0}'.format(pemname)

    # READ THE IMAGE FILE
    #imgname = trainingImgNames[index]
    #intimg  = sitk.ReadImage(imgname)

    # READ THE GENERATED VTK FILES AND FIND THE CORRESPONDING WINDOW FOR PEMS
    vtkname = trainingVtkNames[index]
    vtkpnts = vtk_point_read(vtkname)
    meanpnt = np.mean(vtkpnts,axis=0)
    meanpnt = pemimg.TransformPhysicalPointToIndex([-meanpnt[0],-meanpnt[1],meanpnt[2]])

    # COLLECT THE SAMPLES AROUND THE MEAN POINT
    lower_boun = meanpnt - (downsampleSize * inputSampleSize / 2)
    upper_boun = np.array(pemimg.GetSize()) - (meanpnt + (downsampleSize * inputSampleSize / 2))

    if len([jj for jj in upper_boun if jj<0]) > 0:
        print pemname
        print index
        continue
    if len([jj for jj in lower_boun if jj<0]) > 0:
        print pemname
        print index
        continue

    # CROP THE PEM IMAGE AND USE IT AS INPUT FEATURES
    croppedimg = sitk.DiscreteGaussian(pemimg, smoothWidth)
    croppedimg = sitk.Crop(croppedimg,lower_boun,upper_boun)
    croppedimg = sitk.Shrink(croppedimg, [downsampleSize,downsampleSize,downsampleSize])
    croppedarr = sitk.GetArrayFromImage(croppedimg)

    scipy.misc.imsave('/homes/oo2113/tmp2/{0}.jpg'.format(pemname.split('/')[-1].split('.nii.gz')[0]), croppedarr[inputSampleSize[2]/2,:,:])
    collectedFtrs[index,:numPemFeatures] = croppedarr.flatten()

    # Collect Landmark based features
    combinations   = list(itertools.combinations(range(numLandmarks), 2))
    pair_distances = []
    for cmb in combinations:
        pnt1 = vtkpnts[cmb[0],:]
        pnt2 = vtkpnts[cmb[1],:]
        pair_distances.append(np.sqrt(np.sum(np.power((pnt1 - pnt2), 2))))
    pair_distances = np.array(pair_distances,dtype=np.float64)

    combinations = list(itertools.combinations(range(len(pair_distances)), 2))
    pair_ratios  = []
    for cmb in combinations:
        dist1 = pair_distances[cmb[0]]
        dist2 = pair_distances[cmb[1]]
        pair_ratios.append(dist1/dist2)
    pair_ratios = np.array(pair_ratios,dtype=np.float64)

    lmfeatures = np.concatenate((pair_distances,pair_ratios))
    collectedFtrs[index,numPemFeatures:numPemFeatures+numLMFeatures] = lmfeatures

    # HOG FEATURES
    from skimage.feature import hog
    image = sitk.Crop(pemimg,lower_boun,upper_boun)
    image = sitk.GetArrayFromImage(image)
    hogFtrs = []
    for sliceId in range(inputSampleSize[2] * downsampleSize):
        fd, hog_image = hog(image[sliceId,:,:], orientations=8, pixels_per_cell=(12, 12), cells_per_block=(1, 1), visualise=True)
        hogFtrs.append(fd)
    collectedFtrs[index,numPemFeatures+numLMFeatures:] = np.array(hogFtrs).flatten()


# DEFINE THE REGRESSOR OBJECT AND IT'S PARAMETERS
regressor = RandomForestRegressor(bootstrap=True,
                                  criterion='mse',
                                  max_depth=None,
                                  max_features='auto',
                                  max_leaf_nodes=None,
                                  min_samples_leaf=2,
                                  min_samples_split=2,
                                  min_weight_fraction_leaf=0.0,
                                  n_estimators=16,
                                  n_jobs=-1,
                                  oob_score=False,
                                  random_state=None,
                                  verbose=1,
                                  warm_start=False)

# PERFORM TRAINING THEN INFERENCE ON TESTING DATASET
regressor.fit(collectedFtrs[:numTrainingSamples,:], collectedLabels[:numTrainingSamples])
predictedLabels = regressor.predict(collectedFtrs[numTrainingSamples:,:])
feat_importance = regressor.feature_importances_

print collectedLabels[numTrainingSamples:]
print predictedLabels

print '\n'
imp_pem = np.sum(feat_importance[:numPemFeatures])
imp_LM  = np.sum(feat_importance[numPemFeatures:numPemFeatures+numLMFeatures])
imp_HoG = np.sum(feat_importance[numPemFeatures+numLMFeatures:])
print 'PEM Importance: {0}'.format(imp_pem)
print 'LM Importance: {0}'.format(imp_LM)
print 'HoG Importance: {0}'.format(imp_HoG)


# EVALUATE THE ACCURACY OF PREDICTIONS
mse = np.power(predictedLabels - collectedLabels[numTrainingSamples:],2)
mse = np.sqrt( np.sum(mse[:]) / len(mse.flatten()) )
print 'Total mean square error is {0}'.format(mse)
