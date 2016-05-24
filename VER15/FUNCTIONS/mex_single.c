/*=======================================================================
* mex_single.c - C MEX file used to get the locations of single-center
*     cells from image M (after eliminating multi-center cells)
*------------------------------------------------------------------------
*     usage: [A_add] = mex_single(M, len_lesnum)
*------------------------------------------------------------------------
*     Input arguments:
*       M:      a K x J matrix of elevations
*       n:      number of single-center cells
*
*     Return arguments:
*       A_add:  a array (n x 2) of zero-value (single-center) cells
*------------------------------------------------------------------------
* This is a MEX-file for MATLAB.
* Created by: Phong Le <phongle1@illinois.edu>
* This code may be freely redistributed in its entirety provided that 
  this copyright notice is	not removed.
*========================================================================
*/
#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
  
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *M, *A_add;
  int i, j, p, K, J, n;
  int Outdims[2];               /* dimensions of the output arrays */

  /* Get pointers to inputs, assign default values if not supplied */
  M = (double *)mxGetPr(prhs[0]); 

  /* Get dimensions of input matrix M */
  K = mxGetM(prhs[0]);
  J = mxGetN(prhs[0]);

  n = mxGetScalar(prhs[1]);
  Outdims[0] = 1; Outdims[1] = n*2;   /* Dimension of histo */

  /* Create arrays for the return arguments */
  A_add = mxGetPr(plhs[0] = mxCreateNumericArray(2, Outdims, mxDOUBLE_CLASS, mxREAL));

  p = 0;
  for (j=1; j<J-1; j++){  
    for (i=1 ; i<K-1; i++){
      if (M[K*j+i] != 0) {
        A_add[2*p]  = j+1;
        A_add[2*p+1]= i+1;
        p++;
      }
    }
  }

} /* end mexFunction() */
