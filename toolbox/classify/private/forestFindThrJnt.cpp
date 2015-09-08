/*******************************************************************************
* Piotr's Image&Video Toolbox      Version 3.24
* Copyright 2013 Piotr Dollar.  [pdollar-at-caltech.edu]
* Please email me if you find bugs, or have suggestions or questions!
* Licensed under the Simplified BSD License [see external/bsd.txt]
*******************************************************************************/
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <mex.h>

typedef unsigned int uint32;
#define gini(p) p*p
#define entropy(p) (-p*flog2(float(p)))

// fast approximate log2(x) from Paul Mineiro <paul@mineiro.com>
inline float flog2( float x ) {
  union { float f; uint32_t i; } vx = { x };
  union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
  float y = float(vx.i); y *= 1.1920928955078125e-7f;
  return y - 124.22551499f - 1.498030302f * mx.f
    - 1.72587999f / (0.3520887068f + mx.f);
}

// perform actual computation
void forestFindThrJoint( int H, int N, int F, int L, const float *data,
  const uint32 *hs, const float *ws, const uint32 *order,
  const int split, const float *offsets, const int regSplit,
  uint32 &fid, float &thr, double &gain )
{

  double *Wl, *Wr, *W; float *data1; uint32 *order1;
  int i, j, k, m, j1, j2, h; double vBst, vInit, v, w, wl, wr, g, gl, gr, y_jk, w_y_jk, H_LN, H_RN;
  double imp_reg_left, imp_reg_right, imp_reg_root, imp_cla_left, imp_cla_right, imp_cla_root;
  Wl=new double[H]; Wr=new double[H]; W=new double[H];

  // Regression variable - declaration (malloc) and initialization
  const int nOutputs = L;
  double **sum_left, **sum_right, **sum_root, **mean_left, **mean_right, **mean_root, **var_right, **var_left, **var_root, **sq_sum_left, **sq_sum_right, **sq_sum_root;
  sum_left    = new double* [nOutputs];   sum_right  = new double* [nOutputs];   sum_root  = new double* [nOutputs];
  mean_left   = new double* [nOutputs];   mean_right = new double* [nOutputs];   mean_root = new double* [nOutputs];
  var_left    = new double* [nOutputs];   var_right  = new double* [nOutputs];   var_root  = new double* [nOutputs];
  sq_sum_left = new double* [nOutputs]; sq_sum_right = new double* [nOutputs]; sq_sum_root = new double* [nOutputs];

  for( k=0; k<nOutputs; k++) {
      sum_left[k]    = new double [H]; sum_right[k]    = new double [H]; sum_root[k]    = new double [H];
      mean_left[k]   = new double [H]; mean_right[k]   = new double [H]; mean_root[k]   = new double [H];
      var_left[k]    = new double [H]; var_right[k]    = new double [H]; var_root[k]    = new double [H];
      sq_sum_left[k] = new double [H]; sq_sum_right[k] = new double [H]; sq_sum_root[k] = new double [H];
      for ( m=0; m<H; m++){
        sum_left[k][m]=0;    sum_right[k][m]=0;    sum_root[k][m]=0;
        mean_left[k][m]=0;   mean_right[k][m]=0;   mean_root[k][m]=0;
        var_left[k][m]=0;    var_right[k][m]=0;    var_root[k][m]=0;
        sq_sum_left[k][m]=0; sq_sum_right[k][m]=0; sq_sum_root[k][m]=0;}}

  // initalization
  vBst = vInit = 0; w = 0; fid = 1; thr = 0;
  imp_reg_left = imp_reg_right = imp_reg_root = imp_cla_left = imp_cla_right = imp_cla_root = 0;

  // compute the impurity of the root node both regression and classification
  // classification part
  for( i=0; i<H; i++ ) W[i] = 0;
  for( j=0; j<N; j++ ) { w+=ws[j]; W[hs[j]-1]+=ws[j]; }
  if( split==0 ) { for( i=0; i<H; i++ ) g+=gini(W[i]); imp_cla_root=(1-g/w/w); vBst=vInit=1.0;}
  if( split==1 ) { for( i=0; i<H; i++ ) g+=entropy(W[i]); imp_cla_root=g/w;    vBst=vInit=1.0;}

  // regression part
  for ( k=0; k<nOutputs; k++ ){
    for( j=0; j<N; j++ ){
        h                  = hs[j]-1;
        y_jk               = offsets[j + k*N];
        w_y_jk             = ws[j] * y_jk;
        sum_root[k][h]    += w_y_jk;
        sq_sum_root[k][h] += w_y_jk * y_jk;}
    for ( m=0; m<H; m++){
      mean_root[k][m] = sum_root[k][m] / W[m];
      var_root[k][m]  = sq_sum_root[k][m] / W[m] - mean_root[k][m] * mean_root[k][m];
      imp_reg_root   += (  (W[m]/w) *  (var_root[k][m] / (double)nOutputs) );
    }
  }

  mexPrintf("reg root: %f\n",imp_reg_root);
  mexPrintf("cla root: %f\n",imp_cla_root);
  mexPrintf("0th class prob: %f\n",W[0]/w);
  mexPrintf("1st class prob: %f\n",W[1]/w);
  mexPrintf("total sample: %f\n",w);

  // LOOP OVER FEATURES, THEN THRESHOLD (DATA IS SORTED BY FEATURE VALUE)
  for( i=0; i<F; i++ ) {
    // INITIALIZE THE PARAMETERS AGAIN AFTER EACH FEATURE
    order1=(uint32*) order+i*N; data1=(float*) data+i*size_t(N);
    for ( k=0; k<nOutputs; k++ ){
        for ( m=0; m<H; m++){
          sum_left[k][m]=0;mean_left[k][m]=0;var_left[k][m]=0;sq_sum_left[k][m]=0;
          sum_right[k][m]=sum_root[k][m]; mean_right[k][m]=mean_root[k][m]; var_right[k][m]=var_root[k][m]; sq_sum_right[k][m]=sq_sum_root[k][m];
        }
    }
    for( j=0; j<H; j++ ) { Wl[j]=0; Wr[j]=W[j]; } gl=wl=0; gr=g; wr=w;

    // LOOP OVER THE DATA POINTS
    for( j=0; j<N-1; j++ ) {
      // CLASSIFICATION PART
      j1=order1[j]; j2=order1[j+1]; h=hs[j1]-1;
      if(split==0) {
        // gini = 1-\sum_h p_h^2; v = gini_l*pl + gini_r*pr
        wl+=ws[j1]; gl-=gini(Wl[h]); Wl[h]+=ws[j1]; gl+=gini(Wl[h]);
        wr-=ws[j1]; gr-=gini(Wr[h]); Wr[h]-=ws[j1]; gr+=gini(Wr[h]);
        v = (wl-gl/wl)/w + (wr-gr/wr)/w;
      } else if (split==1) {
        // entropy = -\sum_h p_h log(p_h); v = entropy_l*pl + entropy_r*pr
        gl+=entropy(wl); wl+=ws[j1]; gl-=entropy(wl);
        gr+=entropy(wr); wr-=ws[j1]; gr-=entropy(wr);
        gl-=entropy(Wl[h]); Wl[h]+=ws[j1]; gl+=entropy(Wl[h]);
        gr-=entropy(Wr[h]); Wr[h]-=ws[j1]; gr+=entropy(Wr[h]);
        v = gl/w + gr/w;
      } else if (split==2) {
        // twoing: v = pl*pr*\sum_h(|p_h_left - p_h_right|)^2 [slow if H>>0]
        j1=order1[j]; j2=order1[j+1]; h=hs[j1]-1;
        wl+=ws[j1]; Wl[h]+=ws[j1]; wr-=ws[j1]; Wr[h]-=ws[j1];
        g=0; for( int h1=0; h1<H; h1++ ) g+=fabs(Wl[h1]/wl-Wr[h1]/wr);
        v = - wl/w*wr/w*g*g;
      }
      imp_cla_left=gl/wl; imp_cla_right=gr/wr;
      imp_reg_left=0.0;   imp_reg_right=0.0;

      //REGRESSION PART
      for ( k=0; k<nOutputs; k++){
        y_jk   = offsets[j1 + k*N];
        w_y_jk = ws[j1] * y_jk;
        sum_left[k][h]    += w_y_jk;
        sum_right[k][h]   -= w_y_jk;
        sq_sum_left[k][h] += w_y_jk * y_jk;
        sq_sum_right[k][h]-= w_y_jk * y_jk;}

      for ( k=0; k<nOutputs; k++){
        mean_left[k][h]  = sum_left[k][h] / Wl[h];
        mean_right[k][h] = sum_right[k][h] / Wr[h];
        var_left[k][h]   = sq_sum_left[k][h]  / Wl[h] - mean_left[k][h]  * mean_left[k][h];
        var_right[k][h]  = sq_sum_right[k][h] / Wr[h] - mean_right[k][h] * mean_right[k][h];
        for ( m=0; m<H; m++){
          imp_reg_left  += (  (Wl[m]/wl) *  (var_left[k][m] / (double)nOutputs) );
          imp_reg_right += (  (Wr[m]/wr) *  (var_right[k][m] / (double)nOutputs) );
        }
      }

      // NORMALIZED GAIN COMPUTATION
      H_LN = 0.5 * ((imp_reg_left/imp_reg_root)  + (imp_cla_left/imp_cla_root));
      H_RN = 0.5 * ((imp_reg_right/imp_reg_root) + (imp_cla_right/imp_cla_root));
      v    = H_LN * (wl/w)  + H_RN * (wr/w);

      if( v<vBst && data1[j2]-data1[j1]>=1e-6f ) {
        vBst=v; fid=i+1; thr=0.5f*(data1[j1]+data1[j2]); }

    } // end of threshold search
  } // end of feature search

  // clean up the arrays
  for( k=0; k<nOutputs; k++) {
      delete [] sum_left[k];    delete [] sum_right[k];    delete [] sum_root[k];
      delete [] mean_left[k];   delete [] mean_right[k];   delete [] mean_root[k];
      delete [] var_left[k];    delete [] var_right[k];    delete [] var_root[k];
      delete [] sq_sum_left[k]; delete [] sq_sum_right[k]; delete [] sq_sum_root[k];
  }

  delete [] sum_left;    delete [] sum_right;    delete [] sum_root;
  delete [] mean_left;   delete [] mean_right;   delete [] mean_root;
  delete [] var_left;    delete [] var_right;    delete [] var_root;
  delete [] sq_sum_left; delete [] sq_sum_right; delete [] sq_sum_root;
  delete [] Wl; delete [] Wr; delete [] W; gain = vInit-vBst;

}



// [fid,thr,gain] = mexFunction(data,hs,ws,order,H,split,offsets,regSplit)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int H, N, F, L, split, regSplit; float *data, *ws, thr, *offsets;
  double gain; uint32 *hs, *order, fid;

  data     = (float*) mxGetData(prhs[0]);
  hs       = (uint32*) mxGetData(prhs[1]);
  ws       = (float*) mxGetData(prhs[2]);
  order    = (uint32*) mxGetData(prhs[3]);
  H        = (int) mxGetScalar(prhs[4]);
  split    = (int) mxGetScalar(prhs[5]);
  offsets  = (float*) mxGetData(prhs[6]);
  regSplit = (int) mxGetScalar(prhs[7]);

  N = (int) mxGetM(prhs[0]);
  F = (int) mxGetN(prhs[0]);
  L = (int) mxGetN(prhs[6]);

  forestFindThrJoint(H,N,F,L,data,hs,ws,order,split,offsets,regSplit,fid,thr,gain);
  plhs[0] = mxCreateDoubleScalar(fid);
  plhs[1] = mxCreateDoubleScalar(thr);
  plhs[2] = mxCreateDoubleScalar(gain);
}
