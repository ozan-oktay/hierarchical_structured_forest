/*******************************************************************************
* Licensed under the Simplified BSD License [see external/bsd.txt]
*******************************************************************************/
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <mex.h>

typedef unsigned int uint32;


//affMat = computeAffMat(data2Img(dids1),nImgs);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  uint32 id1, id2, *imgInd;
  int nImgs, nPnts;

  imgInd = (uint32*) mxGetData(prhs[0]);
  nImgs  = (int)   mxGetScalar(prhs[1]);
  nPnts  = (int)   mxGetM(prhs[0]);

  const int outDims[2] = {nImgs,nImgs};

  // create outputs
  plhs[0] = mxCreateNumericArray(2,(const mwSize*)outDims,mxSINGLE_CLASS,mxREAL);
  float *affMat = (float*) mxGetData(plhs[0]);

  // create the matrix
  for (int i=0; i<nPnts; i++){
    for (int j=0; j<nPnts; j++){
        id1 = imgInd[i]-1;
        id2 = imgInd[j]-1;
        affMat[id1+nImgs*id2] += 1.0;
    }
  }
}
