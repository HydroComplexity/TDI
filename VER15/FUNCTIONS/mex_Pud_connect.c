#include "mex.h"
#include <string.h>
#include "matrix.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize M, N, nz, nrow;
  int len_grenum;
  mwIndex i,j,k,h;
  double *pr, *DEM, *linC, *ind_hist;
  mwIndex *ir, *jc;

  if (! mxIsSparse (prhs[0]))
    mexErrMsgTxt ("expects sparse matrix");

  M     = mxGetM (prhs [0]);
  N     = mxGetN (prhs [0]);

  pr    = mxGetPr (prhs[0]);
  ir    = mxGetIr (prhs[0]);
  jc    = mxGetJc (prhs[0]);

  DEM   = (double *)mxGetPr(prhs[1]);
  linC  = (double *)mxGetPr(prhs[2]);

  len_grenum	= mxGetScalar(prhs[3]);
  ind_hist		= (double *)mxGetPr(prhs[4]);

  for (k=1;k<len_grenum;k++){
    h = 0;
    for (j=0;j<N;j++) {
      nrow = jc[j+1] - jc[j];
      for (i = 0; i < nrow; i++){
        if (pr[h] == linC[(int)ind_hist[k]-1])
          DEM[M*j+ (int)ir[h]] = 1;
        h++;
      }
    }
  }
} /* end mexFunction() */
