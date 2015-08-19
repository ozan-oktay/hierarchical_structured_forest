__author__ = 'oo2113'

import os
import irtk
from sklearn.ensemble import RandomForestRegressor
import numpy as np
import SimpleITK as sitk
import vtk
import scipy

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
    points = np.array(points)
    return points

# TRAINING PEM DIRECTORY
trainingImgDir   = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mritestingdatalarge/images'
trainingPemDir   = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mritestingdatalarge/pems_hough_votes_max'
trainingDofDir   = '/vol/biomedic/users/oo2113/str_hier_forest_mri/mritestingdatalarge/dofs'
trainingPemNames = []
trainingDofNames = []
trainingVtkNames = []
trainingImgNames = []

# TRAINING PARAMETERS
numAllSamples      = np.array([900])
numTrainingSamples = np.array([850])
inputSampleSize    = np.array([48,48,8])
downsampleSize     = 2

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
    assert( trainingPemName.split('/')[-1].split('.')[0] == trainingDofName.split('/')[-1].split('.')[0] )
    assert( trainingPemName.split('/')[-1].split('.')[0] == trainingImgName.split('/')[-1].split('.')[0] )
    assert( trainingPemName.split('/')[-1].split('.')[0] == trainingVtkName.split('/')[-1].split('.')[0].split('_lm')[0] )

numAllSamples = np.min([numAllSamples,len(trainingPemNames)])
numFeatures   = 2*inputSampleSize[0]*inputSampleSize[1]*inputSampleSize[2]
collectedLabels = np.zeros([numAllSamples],dtype=float)
collectedPems   = np.zeros([numAllSamples,numFeatures],dtype=float)

# READ THE DOF FILES AND STORE GROUND TRUTH INFORMATION
for index in range(numAllSamples):

    dofname   = trainingDofNames[index]
    dofparams = irtk.AffineTransformation(filename=dofname)
    collectedLabels[index] = dofparams.rz

    # READ THE GENERATED PEM FILE
    pemname = trainingPemNames[index]
    pemimg  = sitk.ReadImage(pemname)

    # READ THE IMAGE FILE
    imgname = trainingImgNames[index]
    intimg  = sitk.ReadImage(imgname)

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
    croppedimg = sitk.Crop(pemimg,lower_boun,upper_boun)
    croppedimg = sitk.DiscreteGaussian(croppedimg, downsampleSize/2)
    croppedimg = sitk.Shrink(croppedimg, [downsampleSize,downsampleSize,downsampleSize])
    croppedarr = sitk.GetArrayFromImage(croppedimg)
    scipy.misc.imsave('/vol/biomedic/users/oo2113/tmp2/{0}.jpg'.format(pemname.split('/')[-1].split('.nii.gz')[0]), croppedarr[4,:,:])
    collectedPems[index,:numFeatures/2] = croppedarr.flatten()

    # CROP THE INTENSITY IMAGE AND USE IT AS ANOTHER SET OF INPUT FEATURES
    croppedint = sitk.Crop(intimg,lower_boun,upper_boun)
    croppedint = sitk.DiscreteGaussian(croppedint, downsampleSize/2)
    croppedint = sitk.Shrink(croppedint, [downsampleSize,downsampleSize,downsampleSize])
    croppedarr = sitk.GetArrayFromImage(croppedint)
    scipy.misc.imsave('/vol/biomedic/users/oo2113/tmp2/{0}_int.jpg'.format(pemname.split('/')[-1].split('.nii.gz')[0]), croppedarr[4,:,:])
    collectedPems[index,numFeatures/2:] = croppedarr.flatten()

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
regressor.fit(collectedPems[:numTrainingSamples,:], collectedLabels[:numTrainingSamples])
predictedLabels = regressor.predict(collectedPems[numTrainingSamples:,:])
feat_importance = regressor.feature_importances_

print collectedLabels[numTrainingSamples:]
print predictedLabels

print '\n'
imp_pem = np.sum(feat_importance[:numFeatures/2])
imp_int = np.sum(feat_importance[numFeatures/2:])
print 'PEM Importance: {0}'.format(imp_pem)
print 'Intensity Importance: {0}'.format(imp_int)

# EVALUATE THE ACCURACY OF PREDICTIONS
mse = np.power(predictedLabels - collectedLabels[numTrainingSamples:],2)
mse = np.sqrt( np.sum(mse[:]) / len(mse.flatten()) )
print 'Total mean square error is {0}'.format(mse)
