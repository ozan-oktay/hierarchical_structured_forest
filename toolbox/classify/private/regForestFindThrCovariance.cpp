#include <string.h>
#include <stdint.h>
#include <math.h>
#include <mex.h>
#define NR_END 1
#define FREE_ARG char*
#define TINY 1.0e-20;
#define TINYIF 1.0e-20

typedef unsigned int uint32;

void nrerror(char const *error_text)
/* Numerical Recipes standard error handler */
{
        fprintf(stderr,"Numerical Recipes run-time error...\n");
        fprintf(stderr,"%s\n",error_text);
        fprintf(stderr,"...now exiting to system...\n");
        exit(1);
}
float *dvector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
        if (!v) nrerror("allocation failure in dvector()");
        return v-nl+NR_END;
}
void free_dvector(float *v, long nl, long nh)
/* free a float vector allocated with dvector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}
int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
        if (!v) nrerror("allocation failure in ivector()");
        return v-nl+NR_END;
}
float **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        float **m;

        /* allocate pointers to rows */
        m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
        if (!m) nrerror("allocation failure 1 in matrix()");
        m += NR_END;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
        if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}
void free_dmatrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by dmatrix() */
{
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}
void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}
void ludcmp(float **a, int n, int *indx, float *d)
/*LU Decomposition - For Determinant Computation */
{
        int i,imax,j,k,m;
        float big,dum,sum,temp;
        float *vv;
        vv=dvector(1,n);
        *d=1.0;
        for (i=1;i<=n;i++) {
                big=0.0;
                for (j=1;j<=n;j++)
                        if ((temp=fabs(a[i][j])) > big) big=temp;
                if (big == 0.0){for (m=1;m<=n;m++) a[m][m]=0.0; return;} //nrerror("Singular matrix in routine ludcmp");
                vv[i]=1.0/big;
        }
        for (j=1;j<=n;j++) {
                for (i=1;i<j;i++) {
                        sum=a[i][j];
                        for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                }
                big=0.0;
                for (i=j;i<=n;i++) {
                        sum=a[i][j];
                        for (k=1;k<j;k++)
                                sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                        if ( (dum=vv[i]*fabs(sum)) >= big) {
                                big=dum;
                                imax=i;
                        }
                }
                if (j != imax) {
                        for (k=1;k<=n;k++) {
                                dum=a[imax][k];
                                a[imax][k]=a[j][k];
                                a[j][k]=dum;
                        }
                        *d = -(*d);
                        vv[imax]=vv[j];
                }
                indx[j]=imax;
                if (a[j][j] == 0.0) a[j][j]=TINY;
                if (j != n) {
                        dum=1.0/(a[j][j]);
                        for (i=j+1;i<=n;i++) a[i][j] *= dum;
                }
        }
        free_dvector(vv,1,n);
}
void Matrix2NR(float **m, float **in, int n)
/*Conversion from float 2D arr to numArray*/
{
  int i, j;
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      m[i+1][j+1] = in[j][i];
    }
  }
}
float Determinant(float **input,int n)
/* LU Decomposition based determinant calculation.*/
{
    float **a,d;
    int i, *index;
    a = dmatrix(1, n, 1, n);
    Matrix2NR(a,input,n);
    index = ivector(1, n);
    ludcmp(a, n, index, &d);
    for (i = 1; i <= n; i++) {
      d *= a[i][i];
    }
    free_dmatrix(a, 1, n, 1, n);
    free_ivector( index, 1, n);
    return d;
}


/*Find the optimal split and return gain,threshold,feature values*/
void forestFindThr( int N, int F, int L, const float *data,
  const float *offsets, const float *ws, const uint32 *order, const int split,
  uint32 &fid, float &thr, float &gain )
{

  // define the parameters - mean and variance parameters
  const int nOutputs = L;
  int i, j, k, m, j1, j2;
  float *data1; uint32 *order1;
  float y_jm, y_jk, w_y_jk, wl, wr, w, vBst, vInit, v, det;
  float impurity_left, impurity_right, impurity_root;

  // IG computation based on determinant of covariance matrix
  float **cov_left, **cov_right, **cov_root, **sq_sum_left, **sq_sum_right, **sq_sum_root;
  float *sum_left, *sum_right, *sum_root, *mean_left, *mean_right, *mean_root;

  cov_left    = new float* [nOutputs]; cov_right    = new float* [nOutputs]; cov_root    = new float* [nOutputs];
  sq_sum_left = new float* [nOutputs]; sq_sum_right = new float* [nOutputs]; sq_sum_root = new float* [nOutputs];
  sum_left    = new float [nOutputs];  sum_right    = new float [nOutputs];  sum_root    = new float [nOutputs];
  mean_left   = new float [nOutputs];  mean_right   = new float [nOutputs];  mean_root   = new float [nOutputs];

  for( m = 0; m < nOutputs; ++m){
         cov_left[m] = new float[nOutputs];    cov_right[m] = new float[nOutputs];    cov_root[m] = new float[nOutputs];
      sq_sum_left[m] = new float[nOutputs]; sq_sum_right[m] = new float[nOutputs]; sq_sum_root[m] = new float[nOutputs];
  }

  // Perform initialization
  vBst = vInit = 0; w = 0; fid = 1; thr = 0;
  impurity_left = impurity_right = impurity_root = 0;
  for( j=0; j<N; j++ ) w+=ws[j];
  for( k=0; k<nOutputs; k++ ){
      sum_left[k]=0;    sum_right[k]=0;    sum_root[k]=0;
      mean_left[k]=0;   mean_right[k]=0;   mean_root[k]=0;
      for ( m=0; m<nOutputs; m++){
          cov_left[k][m]=0;       cov_right[k][m]=0;    cov_root[k][m]=0;
          sq_sum_left[k][m]=0; sq_sum_right[k][m]=0; sq_sum_root[k][m]=0;
      }
  }

  // compute the mean value of the offsets and the cross-product terms
  for ( k=0; k<nOutputs; k++ ){
    for( j=0; j<N; j++ ){
        y_jk            = offsets[j + k*N];
        w_y_jk          = ws[j] * y_jk;
        sum_root[k]    += w_y_jk;
        for ( m=0; m<nOutputs; m++ ){
            y_jm = offsets[j + m*N];
            sq_sum_root[k][m] += w_y_jk * y_jm;
        }
    }
    mean_root[k]   = sum_root[k] / w;
  }

  // compute the impurity of the root node - by computing the covariance matrix
  for ( k=0; k<nOutputs; k++){for ( m=0; m<nOutputs; m++){cov_root[k][m] = sq_sum_root[k][m] / w - mean_root[k] * mean_root[m];}}
  //mexPrintf("root cov:\n",impurity_root);for ( k=0; k<nOutputs; k++){for ( m=0; m<nOutputs; m++){mexPrintf("%f ",cov_root[k][m]);}}mexPrintf("\n");

  det = Determinant (cov_root,nOutputs); if(det<TINYIF) det=HUGE;   impurity_root = log(det);
  //mexPrintf("root det: %f\n",det);

  vBst=vInit=impurity_root;
  //mexPrintf("root impurity: %f\n",impurity_root);

  // loop over features, then thresholds (data is sorted by feature value)
  for( i=0; i<F; i++ ){
    order1=(uint32*) order+i*N; data1=(float*) data+i*size_t(N);
    for ( k=0; k<nOutputs; k++ ) {
        sum_left[k] =0; mean_left[k]=0; sum_right[k]=sum_root[k]; mean_right[k]=mean_root[k];
        for ( m=0; m<nOutputs; m++){
            cov_left[k][m]=0; cov_right[k][m]=0; sq_sum_left[k][m]=0; sq_sum_right[k][m]=sq_sum_root[k][m];
        }
    }
    impurity_left=0; wl=0; impurity_right=impurity_root; wr=w;

    for( j=0; j<N-1; j++ ) {
        j1=order1[j]; j2=order1[j+1];

        for ( k=0; k<nOutputs; k++){
            y_jk   = offsets[j1 + k*N];
            w_y_jk = ws[j1] * y_jk;
            sum_left[k] += w_y_jk;
            sum_right[k]-= w_y_jk;
            for ( m=0; m<nOutputs; m++){
                y_jm = offsets[j1 + m*N];
                sq_sum_left[k][m] += w_y_jk * y_jm;
                sq_sum_right[k][m]-= w_y_jk * y_jm;
            }
        }

        wl+=ws[j1];wr-=ws[j1];v=0;

        for ( k=0; k<nOutputs; k++){ mean_left[k]=sum_left[k]/wl; mean_right[k]=sum_right[k]/wr;}
        for ( k=0; k<nOutputs; k++){
          for ( m=0; m<nOutputs; m++){
              cov_left[k][m]  = sq_sum_left[k][m] / wl  - mean_left[k]  * mean_left[m];
              cov_right[k][m] = sq_sum_right[k][m] / wr - mean_right[k] * mean_right[m];
          }
        }

        det = Determinant (cov_left,nOutputs);  if(det<TINYIF)det=HUGE; impurity_left  = log(det);
        det = Determinant (cov_right,nOutputs); if(det<TINYIF)det=HUGE; impurity_right = log(det);

        v = (wl/w)*impurity_left + (wr/w)*impurity_right;

        if( v<vBst && data1[j2]-data1[j1]>=1e-6f ) {
          vBst=v; fid=i+1; thr=0.5f*(data1[j1]+data1[j2]);}
    }
  }

  // free dynamically allocated memory
  for( i=0; i<nOutputs; i++ ){
      delete[] sq_sum_left[i]; delete[] sq_sum_right[i]; delete[] sq_sum_root[i];// delete array within matrix
      delete[] cov_left[i];    delete[] cov_right[i];    delete[] cov_root[i];// delete array within matrix
  }
  delete [] sum_left;    delete [] sum_right;    delete [] sum_root;
  delete [] mean_left;   delete [] mean_right;   delete [] mean_root;
  delete [] sq_sum_left; delete [] sq_sum_right; delete [] sq_sum_root;
  delete [] cov_left;    delete [] cov_right;    delete [] cov_root;

  gain = vInit-vBst;
  //mexPrintf("best impurity %f\n",vBst);
}


// [fid,thr,gain] = mexFunction(data,offsets,ws,order,split);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int N, F, L, split; float *data, *offsets, *ws, thr;
  float gain; uint32 *order, fid;

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
