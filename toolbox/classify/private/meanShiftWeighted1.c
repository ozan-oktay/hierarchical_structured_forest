/**************************************************************************
 *************************************************************************/
#include <math.h>
#include <string.h>
#include "mex.h"
#ifdef USEOMP
#include <omp.h>
#endif

/**************************************************************************
 * Calculates mean of all the points in data that lie on a sphere of
 * radius^2==radius2 centered on [1xp] vector x. data is [nxp]. mean
 * contains [1xp] result and return is number of points used for calc.
 *************************************************************************/
double meanVec( double *x, double *data, int p, int n, double radius2,
        double *mean, double *weights ) {
  int i, j; double dist; int cnt=0; double m=0;
  for( j=0; j<p; j++ ) mean[j]=0;
  for( i=0; i<n; i++ ) {
    dist = 0.0;
    for( j=0; j<p; j++ ) {
      dist += (x[j]-data[cnt])*(x[j]-data[cnt]); cnt++;
    }
    if( dist < radius2 ) {
      cnt-=p; m+=weights[i];
      for( j=0; j<p; j++ ) mean[j]+=(data[cnt++]*weights[i]);
    }
  }
  if( m ) for( j=0; j<p; j++ ) mean[j]/=m;
  return m;
}

/* Squared euclidean distance between two vectors. */
double dist( double *A, double *B, int n ) {
  double d=0.0; int i;
  for(i=0; i<n; i++) d+=(A[i]-B[i]) * (A[i]-B[i]);
  return d;
}

/**************************************************************************
 * data				- p x n column matrix of data points
 * p                - dimension of data points
 * n                - number of data points
 * radius			- radius of search windo
 * rate				- gradient descent proportionality factor
 * maxIter			- max allowed number of iterations
 * labels			- labels for each cluster
 * means            - output (final clusters)
 *************************************************************************/
void meanShift( double data[], int p, int n, double radius, double rate,
        int maxIter, double labels[], double *means, double *weights, int nThreads ) {
  double radius2;		/* radius^2 */
  int i, j, g;  		/* looping and temporary variables */
  int *deltas;          /* indicator if change occurred between iterations per point */
  double *meansCur;     /* calculated means for current iter */
  double *meansNxt;     /* calculated means for next iter */
  int *consolidated;    /* Needed in the assignment of cluster labels */
  int nLabels = 1;      /* Needed in the assignment of cluster labels */

  /* initialization */
  meansCur = (double*) malloc( sizeof(double)*p*n );
  meansNxt = (double*) malloc( sizeof(double)*p*n );
  consolidated = (int*) malloc( sizeof(int)*n );
  deltas = (int*) malloc( sizeof(int)*n );
  for(i=0; i<n; i++) deltas[i] = 1;
  radius2 = radius * radius;
  meansCur = (double*) memcpy(meansCur, data, p*n*sizeof(double) );

  /* main loop */
  #ifdef USEOMP
  nThreads = (nThreads < omp_get_max_threads()) ? nThreads : omp_get_max_threads();
  #pragma omp parallel for num_threads(nThreads) shared(data,weights,meansCur,meansNxt) private(i)
  #endif
  for ( i=0; i<n; i++) {
    int k,iter,o=i*p; double m,*mean;
    mean = (double*) malloc( sizeof(double)*p );
    for( iter=0; iter<maxIter; iter++ ) {
      if( deltas[i] ) {
        /* shift meansNxt in direction of mean (if m>0) */
        m=meanVec( meansCur+o, data, p, n, radius2, mean, weights );
        if( m ) {
          for( k=0; k<p; k++ ) meansNxt[o+k] = (1-rate)*meansCur[o+k] + rate*mean[k];
          if( dist(meansNxt+o, meansCur+o, p)<0.00003) deltas[i]=0;
          for( k=0; k<p; k++ ) meansCur[o+k] = meansNxt[o+k];
        } else {
          deltas[i]=0;
        }
      }
    }
    free(mean);
  }

  /* Consolidate: assign all points that are within radius2 to same cluster. */
  for( i=0; i<n; i++ ) { consolidated[i]=0; labels[i]=0; }
  for( i=0; i<n; i++ ) if( !consolidated[i]) {
    for( j=0; j<n; j++ ) if( !consolidated[j]) {
      if( dist(meansCur+i*p, meansCur+j*p, p) < radius2) {
        labels[j]=nLabels; consolidated[j]=1;
      }
    }
    nLabels++;
  }
  nLabels--; memcpy( means, meansCur, p*n*sizeof(double) );

  /* free memory */
  free(meansNxt); free(meansCur); free(consolidated); free(deltas);
}

/* see meanShift.m for usage info */
void mexFunction( int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[] ) {
  double radius, rate, *data, *labels, *means, *weights; int p, n, maxIter, nThreads;

  /* Check inputs */
  if(nrhs < 4) mexErrMsgTxt("At least four input arguments required.");
  if(nlhs > 2) mexErrMsgTxt("Too many output arguments.");

  /* Get inputs */
  data = mxGetPr(prhs[0]);
  radius = mxGetScalar(prhs[1]);
  rate = mxGetScalar(prhs[2]);
  maxIter = (int) mxGetScalar(prhs[3]);
  weights = mxGetPr(prhs[4]);
  nThreads = (int) mxGetScalar(prhs[5]);
  p=mxGetM(prhs[0]); n=mxGetN(prhs[0]);

  /* Create outputs */
  plhs[0] = mxCreateNumericMatrix(n, 1, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericMatrix(p, n, mxDOUBLE_CLASS, mxREAL);
  labels=mxGetPr(plhs[0]); means=mxGetPr(plhs[1]);

  /* Do the actual computations in a subroutine */
  meanShift( data, p, n, radius, rate, maxIter, labels, means, weights, nThreads );
}