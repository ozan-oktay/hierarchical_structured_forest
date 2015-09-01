/*******************************************************************************
* Code written by Ozan Oktay, 2015.
* Licensed under the MSR-LA Full Rights License [see license.txt]
*******************************************************************************/
#include <mex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef USEOMP
#include <omp.h>
#endif

typedef unsigned long long int 	uint64;
typedef unsigned int uint32;
typedef unsigned short uint16;
typedef unsigned char uint8;
template<typename T> inline T min( T x, T y ) { return x<y ? x : y; }

int round_int( float r ) {
    return (r > 0.0) ? (r + 0.5) : (r - 0.5);
}

// construct lookup array for mapping fids to channel indices
uint32* buildLookup( int *dims, int w ) {
  int c, r, s, z, n=w*w*w*dims[3]; uint32 *cids=new uint32[n]; n=0;
  for(z=0; z<dims[3]; z++)
      for(s=0; s<w; s++)
          for(c=0; c<w; c++)
              for(r=0; r<w; r++)
                  cids[n++] = z*dims[0]*dims[1]*dims[2]
                            + s*dims[0]*dims[1]
                            + c*dims[0]
                            + r;
  return cids;
}

// construct lookup arrays for mapping fids for self-similarity channel
void buildLookupSs( uint32 *&cids1, uint32 *&cids2, int *dims, int w, int m ) {
  int i, j, z, z1, d, c, r; int locs[1024];
  int m3=m*m*m, n=m3*(m3-1)/2*dims[3], s=int(w/m/2.0+.5);
  cids1 = new uint32[n]; cids2 = new uint32[n]; n=0;
  for(i=0; i<m; i++) locs[i]=uint32((i+1)*(w+2*s-1)/(m+1.0)-s+.5);
  for(z=0; z<dims[3]; z++)
      for(i=0; i<m3; i++)
          for(j=i+1; j<m3; j++) {
              z1=z*dims[0]*dims[1]*dims[2]; n++;
              r=i%m; c=((i-r)/m)%m; d=(i-(r+c*m))/(m*m); cids1[n-1]= z1 + locs[d]*dims[0]*dims[1] + locs[c]*dims[0] + locs[r];
              r=j%m; c=((j-r)/m)%m; d=(j-(r+c*m))/(m*m); cids2[n-1]= z1 + locs[d]*dims[0]*dims[1] + locs[c]*dims[0] + locs[r];
          }
}

// [E,ind,HV] = mexFunction(model,I,chns,chnsSs,chnsSh) - helper for edgesDetect.m
void mexFunction( int nl, mxArray *pl[], int nr, const mxArray *pr[] )
{
  // get inputs
  mxArray *model = (mxArray*) pr[0];
  float *I = (float*) mxGetData(pr[1]);
  float *chns = (float*) mxGetData(pr[2]);
  float *chnsSs = (float*) mxGetData(pr[3]);
  float *chnsSh = (float*) mxGetData(pr[4]);
  bool chnsShEmpty = mxIsEmpty(pr[4]);

  // extract relevant fields from model and options
  float *thrs = (float*) mxGetData(mxGetField(model,0,"thrs"));
  uint32 *fids = (uint32*) mxGetData(mxGetField(model,0,"fids"));
  uint32 *child = (uint32*) mxGetData(mxGetField(model,0,"child"));
  uint8 *segs = (uint8*) mxGetData(mxGetField(pr[0],0,"segs"));
  uint8 *nSegs = (uint8*) mxGetData(mxGetField(pr[0],0,"nSegs"));
  uint16 *eBins = (uint16*) mxGetData(mxGetField(model,0,"eBins"));
  uint32 *eBnds = (uint32*) mxGetData(mxGetField(model,0,"eBnds"));
  mxArray *opts = mxGetField(model,0,"opts");

  float *meanOff = (float*) mxGetData(mxGetField(model,0,"meanOff"));
  float *covOff = (float*) mxGetData(mxGetField(model,0,"covOff"));
  float *meanPose = (float*) mxGetData(mxGetField(model,0,"meanPose"));
  float *covPose = (float*) mxGetData(mxGetField(model,0,"covPose"));

  const int shrink = (int) mxGetScalar(mxGetField(opts,0,"shrink"));
  const int imWidth = (int) mxGetScalar(mxGetField(opts,0,"imWidth"));
  const int gtWidth = (int) mxGetScalar(mxGetField(opts,0,"gtWidth"));
  const int nChns = (int) mxGetScalar(mxGetField(opts,0,"nChns"));
  const int nCells = (int) mxGetScalar(mxGetField(opts,0,"nCells"));
  const int nLandmarks = (int) mxGetScalar(mxGetField(opts,0,"nLandmarks"));
  const int stageId = (int) mxGetScalar(mxGetField(opts,0,"stageId"));
  const int nPosePar = (int) mxGetScalar(mxGetField(opts,0,"nPosePar"));
  const uint32 nChnFtrs = (uint32) mxGetScalar(mxGetField(opts,0,"nChnFtrs"));
  const uint32 nSimFtrs = (uint32) mxGetScalar(mxGetField(opts,0,"nSimFtrs"));
  const int stride = (int) mxGetScalar(mxGetField(opts,0,"stride"));
  const int nTreesEval = (int) mxGetScalar(mxGetField(opts,0,"nTreesEval"));
  int sharpen = (int) mxGetScalar(mxGetField(opts,0,"sharpen"));
  int nThreads = (int) mxGetScalar(mxGetField(opts,0,"nThreads"));
  const int nBnds = int(mxGetNumberOfElements(mxGetField(model,0,"eBnds"))-1)/
    int(mxGetNumberOfElements(mxGetField(model,0,"thrs")));
  const char *msgSharpen="Model supports sharpening of at most %i pixels!\n";
  if( sharpen>nBnds-1 ) { sharpen=nBnds-1; mexPrintf(msgSharpen,sharpen); }
  if ( (stageId == 1) & chnsShEmpty) mexErrMsgTxt("Channel Features shouldn't be empty for hier level=1.");

  // get dimensions and constants
  const mwSize *imgSize = mxGetDimensions(pr[1]);
  const int h = (int) imgSize[0];
  const int w = (int) imgSize[1];
  const int d = (int) imgSize[2];
  const int Z = mxGetNumberOfDimensions(pr[1])<=3 ? 1 : imgSize[3];
  const mwSize *fidsSize = mxGetDimensions(mxGetField(model,0,"fids"));
  const int nTreeNodes = (int) fidsSize[0];
  const int nTrees = (int) fidsSize[1];
  const int h1 = (int) ceil(double(h-imWidth)/stride);
  const int w1 = (int) ceil(double(w-imWidth)/stride);
  const int d1 = (int) ceil(double(d-imWidth)/stride);
  const int h2 = h1*stride+gtWidth;
  const int w2 = w1*stride+gtWidth;
  const int d2 = d1*stride+gtWidth;
  const int hCov = nLandmarks * 3;
  const int hOff = nLandmarks * 3;
  const int imgDims[4] = {h,w,d,Z};
  const int chnDims[4] = {h/shrink,w/shrink,d/shrink,nChns};
  const int indDims[4] = {h1,w1,d1,nTreesEval};
  const int outDims[4] = {h2,w2,d2,1};
  const int hovDims[4] = {h2,w2,d2,nLandmarks};
  const int posDims[2] = {2,nPosePar};
  const int imRadius   = imWidth/2;
  const int gtRadius   = gtWidth/2;

  // construct lookup tables
  uint32 *iids, *eids, *cids, *cids1, *cids2;
  iids = buildLookup( (int*)imgDims, gtWidth );
  eids = buildLookup( (int*)outDims, gtWidth );
  cids = buildLookup( (int*)chnDims, imWidth/shrink );
  buildLookupSs( cids1, cids2, (int*)chnDims, imWidth/shrink, nCells );

  // create outputs
  pl[0] = mxCreateNumericArray(4,(const mwSize*)outDims,mxSINGLE_CLASS,mxREAL);
  float *E = (float*) mxGetData(pl[0]);
  pl[1] = mxCreateNumericArray(4,(const mwSize*)indDims,mxUINT32_CLASS,mxREAL);
  uint32 *ind = (uint32*) mxGetData(pl[1]);
  pl[2] = mxCreateNumericArray(4,(const mwSize*)hovDims,mxSINGLE_CLASS,mxREAL);
  float *HV = (float*) mxGetData(pl[2]);
  pl[3] = mxCreateNumericArray(2,(const mwSize*)posDims,mxSINGLE_CLASS,mxREAL);
  float *P = (float*) mxGetData(pl[3]);

  // apply forest to all patches and store leaf inds
  #ifdef USEOMP
  nThreads = min(nThreads,omp_get_max_threads());
  #pragma omp parallel for num_threads(nThreads)
  #endif
  for (int sl=0; sl<d1; sl++) for( int c=0; c<w1; c++ ) for( int t=0; t<nTreesEval; t++ ) {
    for( int r0=0; r0<2; r0++ ) for( int r=r0; r<h1; r+=2 ) {
      int o = (r*stride/shrink) + (c*stride/shrink)*h/shrink + (sl*stride/shrink)*(h/shrink)*(w/shrink);
      // select tree to evaluate
      int t1 = ((r+c+sl)%2*nTreesEval+t)%nTrees; uint32 k = t1*nTreeNodes;
      while( child[k] ) {
        // compute feature (either channel or self-similarity feature)
        uint32 f = fids[k]; float ftr;
        if( f<nChnFtrs ) ftr = chns[cids[f]+o]; // channel features
        else if ( f<(nChnFtrs+nSimFtrs) ) ftr = chnsSs[cids1[f-nChnFtrs]+o]-chnsSs[cids2[f-nChnFtrs]+o]; // self-similarity features
        else ftr = chnsSh[f-(nChnFtrs+nSimFtrs)]; // shape features
        // compare ftr to threshold and move left or right accordingly
        if( ftr < thrs[k] ) k = child[k]-1; else k = child[k];
        k += t1*nTreeNodes;
      }
      // store leaf index and update edge maps
      ind[ r + c*h1 + sl*h1*w1 + t*h1*w1*d1 ] = k;
    }
  }

  // compute edge maps (avoiding collisions from parallel executions)
  if( !sharpen ) for( int c0=0; c0<gtWidth/stride; c0++ ) {
  #ifdef USEOMP
  #pragma omp parallel for num_threads(nThreads)
  #endif
    for( int c=c0; c<w1; c+=gtWidth/stride ) {
      for( int t=0; t<nTreesEval; t++ ) {
        for ( int sl=0; sl<d1; sl++ ) {
          for( int r=0; r<h1; r++ ) {
            uint32 k = ind[ r + c*h1 + sl*h1*w1 + t*h1*w1*d1 ];
            float *E1 = E + (r*stride) + (c*stride)*h2 + (sl*stride)*h2*w2;
            int b0=eBnds[k*nBnds], b1=eBnds[k*nBnds+1]; if(b0==b1) continue;
            for( int b=b0; b<b1; b++ ) E1[eids[eBins[b]]]++;
          }
        }
      }
    }
  }

  // computed sharpened edge maps, snapping to local color values
  if( sharpen ) {
    // compute neighbors array
    const int g=gtWidth; uint32 N[65536*6];
    for( int s=0; s<g; s++ ) for( int c=0; c<g; c++ ) for( int r=0; r<g; r++ ) {
      int i=s*g*g+c*g+r; uint32 *N1=N+i*6;
      N1[0] = c>0 ? i-g   : i; N1[1] = c<g-1 ? i+g   : i;
      N1[2] = r>0 ? i-1   : i; N1[3] = r<g-1 ? i+1   : i;
      N1[4] = s>0 ? i-g*g : i; N1[5] = s<g-1 ? i+g*g : i;
    }
    #ifdef USEOMP
    #pragma omp parallel for num_threads(nThreads)
    #endif
    for( int c=0; c<w1; c++ ) for ( int sl=0; sl<d1; sl++ ) for( int r=0; r<h1; r++ ) {
      for( int t=0; t<nTreesEval; t++ ) {
        // get current segment and copy into S
        uint64 k   = ind[ r + c*h1 + sl*h1*w1 + t*h1*w1*d1 ];
        uint64 g_l = g;
        int m = nSegs[k]; if( m==1 ) continue;
        uint8 S0[65536], *S=S0;
        memcpy(S,segs+k*g_l*g_l*g_l, g*g*g*sizeof(uint8));
        // compute color model for each segment using every other pixel
        int si, ci, ri, s, z; float ns[100], mus[1000];
        const float *I1 = I+(sl*stride+(imWidth-g)/2)*h*w+(c*stride+(imWidth-g)/2)*h+r*stride+(imWidth-g)/2;
        for( s=0; s<m; s++ ) { ns[s]=0; for( z=0; z<Z; z++ ) mus[s*Z+z]=0; }
        for( si=0; si<g; si+=2 ) for( ci=0; ci<g; ci+=2 ) for( ri=0; ri<g; ri+=2 ) {
          s = S[si*g*g+ci*g+ri]; ns[s]++;
          for( z=0; z<Z; z++ ) mus[s*Z+z]+=I1[z*h*w*d+si*h*w+ci*h+ri];
        }
        for(s=0; s<m; s++) for( z=0; z<Z; z++ ) mus[s*Z+z]/=ns[s];

        // update segment S according to local color values
        int b0=eBnds[k*nBnds], b1=eBnds[k*nBnds+sharpen];
        for( int b=b0; b<b1; b++ ) {
          float vs[10], dx, e, eBest=1e10f; int i, sBest=-1, ss[6];
          for( i=0; i<6; i++ ) ss[i]=S[N[eBins[b]*6+i]];
          for( z=0; z<Z; z++ ) vs[z]=I1[iids[eBins[b]]+z*h*w*d];
          for( i=0; i<6; i++ ) {
            s=ss[i]; if(s==sBest) continue;
            e=0; for( z=0; z<Z; z++ ) { dx=mus[s*Z+z]-vs[z]; e+=dx*dx; }
            if( e<eBest ) { eBest=e; sBest=s; }
          }
          S[eBins[b]]=sBest;
        }
        // convert mask to edge maps (examining expanded set of pixels)
        float *E1 = E + sl*stride*h2*w2 + c*stride*h2 + r*stride; b1=eBnds[k*nBnds+sharpen+1];
        for( int b=b0; b<b1; b++ ) {
          int i=eBins[b]; uint8 s=S[i]; uint32 *N1=N+i*6;
          if( s!=S[N1[0]] || s!=S[N1[1]] || s!=S[N1[2]] || s!=S[N1[3]] || s!=S[N1[4]] || s!=S[N1[5]] )
            E1[eids[i]]++;
        }
      }
    }
  }

  // compute the pose parameters
  if (stageId==1 && nPosePar>0){
    for( int t=0; t<nTreesEval; t++ ) {
        for ( int sl=0; sl<d1; sl++ ) {
            for( int c=0; c<w1; c++ ) {
                for( int r=0; r<h1; r++ ) {
                    uint32 k = ind[ r + c*h1 + sl*h1*w1 + t*h1*w1*d1 ];
                    for ( int p=0; p<nPosePar; p++){
                        P[p*2]   += meanPose[k*nPosePar + p];
                        P[p*2+1] += covPose [k*nPosePar*nPosePar + p*nPosePar + p];
                    }
                }
            }
        }
    }
    for ( int p=0; p<nPosePar*2; p++) P[p] /= (float)(nTreesEval*d1*w1*h1);
  }

  // compute the mean confidence value
  uint32 vCount=0; float *covMean; covMean = new float [nLandmarks]; for (int lm=0; lm<nLandmarks; lm++) covMean[lm]=0.0f;
  if (nLandmarks > 0){
      for( int t=0; t<nTreesEval; t++ ) {
          for ( int sl=0; sl<d1; sl++ ) {
              for( int c=0; c<w1; c++ ) {
                  for( int r=0; r<h1; r++ ) {
                      uint32 k = ind[ r + c*h1 + sl*h1*w1 + t*h1*w1*d1 ];
                      for ( int lm=0; lm<nLandmarks; lm++) {
                        float covR = covOff[ (lm*3)  *hCov + (lm*3)   + k*hCov*hCov ];
                        float covC = covOff[ (lm*3+1)*hCov + (lm*3+1) + k*hCov*hCov ];
                        float covS = covOff[ (lm*3+2)*hCov + (lm*3+2) + k*hCov*hCov ];
                        float covTr = (covR + covC + covS);
                        covMean[lm] += covTr; if(covTr>1e-15) vCount++;
                      }
                  }
              }
          }
      }
      for (int lm=0; lm<nLandmarks; lm++) covMean[lm] /= ( (float)vCount / (float)nLandmarks );
      for (int lm=0; lm<nLandmarks; lm++) mexPrintf("covMean %d: %f\n",lm,covMean[lm]);
  }

  // compute Hough Votes for the landmarks
  if (nLandmarks > 0){
      #ifdef USEOMP
      #pragma omp parallel for num_threads(nThreads)
      #endif
      for( int t=0; t<nTreesEval; t++ ) {
          for ( int sl=0; sl<d1; sl++ ) {
              for( int c=0; c<w1; c++ ) {
                  for( int r=0; r<h1; r++ ) {

                      uint32 k = ind[ r + c*h1 + sl*h1*w1 + t*h1*w1*d1 ];

                      int ri  = r*stride  + gtRadius - 1;
                      int ci  = c*stride  + gtRadius - 1;
                      int si  = sl*stride + gtRadius - 1;
                      float E1  = *(E + ri + ci*h2 + si*h2*w2);

                      for ( int lm=0; lm<nLandmarks; lm++) {

                        float ro = meanOff[     lm*3 + k*hOff ];
                        float co = meanOff[ 1 + lm*3 + k*hOff ];
                        float so = meanOff[ 2 + lm*3 + k*hOff ];

                        float covR = covOff[ (lm*3)  *hCov + (lm*3)   + k*hCov*hCov ];
                        float covC = covOff[ (lm*3+1)*hCov + (lm*3+1) + k*hCov*hCov ];
                        float covS = covOff[ (lm*3+2)*hCov + (lm*3+2) + k*hCov*hCov ];
                        float covTr = (covR + covC + covS);
                        if (covTr >= covMean[lm]/4.0f) continue;

                        float rn = (float)ri + ro; float rn_p = rn - (float)floor(rn); int rn_f = int(floor(rn)); int rn_c = int(ceil(rn));
                        float cn = (float)ci + co; float cn_p = cn - (float)floor(cn); int cn_f = int(floor(cn)); int cn_c = int(ceil(cn));
                        float sn = (float)si + so; float sn_p = sn - (float)floor(sn); int sn_f = int(floor(sn)); int sn_c = int(ceil(sn));

                        if (rn_f<0 || cn_f<0 || sn_f<0 || rn_c>=h2 || cn_c>=w2 || sn_c>=d2) continue;

                        uint64 hv1 = rn_f + cn_f*h2 + sn_f*h2*w2 + lm*h2*w2*d2; HV[hv1] += E1*(1-rn_p)*(1-cn_p)*(1-sn_p);
                        uint64 hv2 = rn_c + cn_f*h2 + sn_f*h2*w2 + lm*h2*w2*d2; HV[hv2] += E1*  (rn_p)*(1-cn_p)*(1-sn_p);
                        uint64 hv3 = rn_f + cn_c*h2 + sn_f*h2*w2 + lm*h2*w2*d2; HV[hv3] += E1*(1-rn_p)*  (cn_p)*(1-sn_p);
                        uint64 hv4 = rn_c + cn_c*h2 + sn_f*h2*w2 + lm*h2*w2*d2; HV[hv4] += E1*  (rn_p)*  (cn_p)*(1-sn_p);

                        uint64 hv5 = rn_f + cn_f*h2 + sn_c*h2*w2 + lm*h2*w2*d2; HV[hv5] += E1*(1-rn_p)*(1-cn_p)*  (sn_p);
                        uint64 hv6 = rn_c + cn_f*h2 + sn_c*h2*w2 + lm*h2*w2*d2; HV[hv6] += E1*  (rn_p)*(1-cn_p)*  (sn_p);
                        uint64 hv7 = rn_f + cn_c*h2 + sn_c*h2*w2 + lm*h2*w2*d2; HV[hv7] += E1*(1-rn_p)*  (cn_p)*  (sn_p);
                        uint64 hv8 = rn_c + cn_c*h2 + sn_c*h2*w2 + lm*h2*w2*d2; HV[hv8] += E1*  (rn_p)*  (cn_p)*  (sn_p);

                      }
                  }
              }
          }
      }
  }

  // free memory
  delete [] iids; delete [] eids;
  delete [] cids; delete [] cids1; delete [] cids2;
  delete [] covMean;
}
