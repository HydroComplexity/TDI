/*=======================================================================
* mex_fill_pud.c - C MEX file used to fill all depression cells to 
      an user-defined value
*------------------------------------------------------------------------
*     usage: [DEM_pud] = mex_fill_pud(DEM, Pud_pts, ind_pudcum)

*------------------------------------------------------------------------
*     Input arguments:
*       DEM:        a M x N matrix of elevation
*       Pud_pts:    List of depression cells
*       ind_pudcum: Cumulative index of ind_pud
*
*     Return arguments:
*       DEM_pud:    a new M x N matrix of elevation
*------------------------------------------------------------------------
* This is a MEX-file for MATLAB.
* Created by: Phong Le <phongle1@illinois.edu>
* This code may be freely redistributed in its entirety provided that 
  this copyright notice is	not removed.
*========================================================================
*/

#include "mex.h"
#include <string.h>
#include "matrix.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *DEM, *Psi, *Pud_pts, *ind_pudcum;
  int i, j, M, N, Mp, num_pud;

  /* Get value and dimensions of input matrix M */
  DEM         = (double *)mxGetPr(prhs[0]);
  M           = mxGetM(prhs[0]);
  N           = mxGetN(prhs[0]);

  Pud_pts     = (double *)mxGetPr(prhs[1]);
  Mp          = mxGetM(prhs[1]);
  ind_pudcum  = (double *)mxGetPr(prhs[2]);
  num_pud     = mxGetNumberOfElements(prhs[2]);

  Psi         = (double*)malloc(M*N*sizeof(double));
  memcpy(Psi, DEM, M*N*sizeof(double));	

  for (i=1;i<num_pud;i++){
    for (j=(int)ind_pudcum[i-1];j<(int)ind_pudcum[i];j++){
      Psi[M*(int)Pud_pts[0*Mp+j]+(int)Pud_pts[1*Mp+j]] = 10000;
    }
  }

  plhs[0] = mxCreateNumericMatrix(M,N, mxDOUBLE_CLASS, mxREAL);
  memmove(mxGetData(plhs[0]),Psi,M*N*sizeof(double));

} /* end mexFunction() */
