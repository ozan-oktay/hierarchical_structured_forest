/*******************************************************************************
* Author : Ozan Oktay
* Date   : July 2015
*******************************************************************************/
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <mex.h>

typedef unsigned int uint32;

// perform actual computation
void forestFindThr( int N, int F, int L, const float *data,
  const float *offsets, const float *ws, const uint32 *order, const int split,
  uint32 &fid, float &thr, double &gain )
{

  // define the parameters - mean and variance parameters
  const int nOutputs = L;
  int i, j, k, j1, j2;
  float *data1; uint32 *order1;
  double y_ik, w_y_ik, w_yy_ik, wl, wr, w, vBst, vInit, v, vl, vr;
  double *sum_left, *sum_right, *sum_root, *sq_sum_left, *sq_sum_right, *sq_sum_root;
  double *impurity_left, *impurity_right, *impurity_root;

  sum_left      = new double [nOutputs];   sum_right    = new double [nOutputs];   sum_root    = new double [nOutputs];
  sq_sum_left   = new double [nOutputs]; sq_sum_right   = new double [nOutputs]; sq_sum_root   = new double [nOutputs];
  impurity_left = new double [nOutputs]; impurity_right = new double [nOutputs]; impurity_root = new double [nOutputs];

  // perform initialization
  vBst = vInit = 0; w = 0; fid = 1; thr = 0;
  for( j=0; j<N; j++ ) w+=ws[j];
  for( k=0; k<nOutputs; k++ ){
      sum_left[k]=0;      sum_right[k]=0;      sum_root[k]=0;
      sq_sum_left[k]=0;   sq_sum_right[k]=0;   sq_sum_root[k]=0;
      impurity_root[k]=0; impurity_right[k]=0; impurity_left[k]=0;
  }

  // compute the impurity of the root node
  for ( k=0; k<nOutputs; k++ ){
    for( j=0; j<N; j++ ){
        y_ik            = offsets[j + k*N];
        w_y_ik          = ws[j] * y_ik;
        w_yy_ik         = w_y_ik * y_ik;
        sum_root[k]    += w_y_ik;
        sq_sum_root[k] += w_yy_ik;
    }
    impurity_root[k] = (sq_sum_root[k] - sum_root[k] * sum_root[k] / w) / w;
  }
  vBst=vInit=1.0;

  // loop over features, then thresholds (data is sorted by feature value)
  for( i=0; i<F; i++ ){
    order1=(uint32*) order+i*N; data1=(float*) data+i*size_t(N);
    for ( k=0; k<nOutputs; k++ ) {
        sum_left[k] =0;            sq_sum_left[k]=0;
        sum_right[k]=sum_root[k];  sq_sum_right[k]=sq_sum_root[k];
    }
    wl=0; wr=w;

    for( j=0; j<N-1; j++ ) {
        j1=order1[j]; j2=order1[j+1];

        for ( k=0; k<nOutputs; k++){
            y_ik    = offsets[j1 + k*N];
            w_y_ik  = ws[j1] * y_ik;
            w_yy_ik = w_y_ik * y_ik;
            sum_left[k] += w_y_ik;
            sum_right[k]-= w_y_ik;
            sq_sum_left[k] += w_yy_ik;
            sq_sum_right[k]-= w_yy_ik;
        }

        wl+=ws[j1];wr-=ws[j1];v=0;vl=0;vr=0;

        for ( k=0; k<nOutputs; k++){
            impurity_left[k]  = (sq_sum_left[k]  - sum_left[k]  * sum_left[k]  / wl) / wl;
            impurity_right[k] = (sq_sum_right[k] - sum_right[k] * sum_right[k] / wr) / wr;
            vl += (impurity_left[k] /impurity_root[k]);
            vr += (impurity_right[k]/impurity_root[k]);
        }
        vl /= (double)nOutputs;
        vr /= (double)nOutputs;
        v   = (wl/w)*vl + (wr/w)*vr;

        if( v<vBst && data1[j2]-data1[j1]>=1e-6f ) {
          vBst=v; fid=i+1; thr=0.5f*(data1[j1]+data1[j2]);}
    }
  }

  // clean up the arrays
  delete [] sum_left;      delete [] sum_right;      delete [] sum_root;
  delete [] sq_sum_left;   delete [] sq_sum_right;   delete [] sq_sum_root;
  delete [] impurity_left; delete [] impurity_right; delete [] impurity_root;
  gain = vInit-vBst;

}


// [fid,thr,gain] = mexFunction(data,offsets,ws,order,split);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int N, F, L, split; float *data, *offsets, *ws, thr;
  double gain; uint32 *order, fid;

  data    = (float*) mxGetData(prhs[0]);
  offsets = (float*) mxGetData(prhs[1]);
  ws      = (float*) mxGetData(prhs[2]);
  order   = (uint32*) mxGetData(prhs[3]);
  split   = (int) mxGetScalar(prhs[4]);

  N = (int) mxGetM(prhs[0]);
  F = (int) mxGetN(prhs[0]);
  L = (int) mxGetN(prhs[1]);

  forestFindThr(N,F,L,data,offsets,ws,order,split,fid,thr,gain);

  plhs[0] = mxCreateDoubleScalar(fid);
  plhs[1] = mxCreateDoubleScalar(thr);
  plhs[2] = mxCreateDoubleScalar(gain);
}
