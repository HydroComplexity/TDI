/*=======================================================================
* mex_fill_sing.c - C MEX file used to fill single-center cells to  value
*     the maximum elevation
*------------------------------------------------------------------------
*     usage: [A] = mex_fill_sing(Center, M, maxE)
*------------------------------------------------------------------------
*     Input arguments:
*       Center: list (x&y) locations of single-center cells
*       M:      a K x J matrix of elevations
*       maxE:   maximum elevation value
*
*     Return arguments:
*       A:      a matrix of new elevations
*------------------------------------------------------------------------
* This is a MEX-file for MATLAB.
* Created by: Phong Le <phongle1@illinois.edu>
*========================================================================
*/
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize n;
  double *M, *A;
  int i, j, K, J, row, col;
  double *vri, *max;

  /* Get pointers to inputs, assign default values if not supplied */
  M = (double *)mxGetPr(prhs[1]); 

  /* Get dimensions of input matrix M */
  n   = mxGetNumberOfElements(prhs[0]);
  vri = mxGetPr(prhs[0]);

  K   = mxGetM(prhs[1]);
  J   = mxGetN(prhs[1]);

  max = mxGetPr(prhs[2]);

  /* Create arrays for the return arguments */
  A = (double *)mxGetPr(plhs[0] = mxDuplicateArray(prhs[1])); /* Duplicate matrix */

  /* Update max elevation to matrix M */
  for (i=0; i<n/2; i++) {
    col = (int) vri[i];
    row = (int) vri[i+n/2];
    A[(col-1)*K + (row-1)] = max[0];
  }
} /* end mexFunction() */
