/*=======================================================================
* mex_geometric.c - C MEX file used to calculate Volume and Depth of     
      depression cells
*------------------------------------------------------------------------
* usage: [volume,depth] = 
      mex_geometric(DEM,Pud_pts,Thres_pts,ind_pudcum,ind_threscum,dx,dy);

*------------------------------------------------------------------------
*     Input arguments:
*       DEM:          a K x J matrix of elevation
*       Pud_pts:      List of depression cells
*       Thres_pts:    List of depression cells
*       ind_pudcum:   Cumulative index of ind_pud
*       ind_threscum: Cumulative index of ind_thres
*
*     Return arguments:
*       Volume:       list of volume of depressions
*       Depth:        list of depth of depressions
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
  double *DEM, *Center_real, *Pud_pts, *Thres_pts, *ind_pudcum, *ind_threscum, *volume, *depth;
  double dx, dy;
  int i, j, k, M, N, Mp, Mt, num_thres, num_pud;

  /* Get value and dimensions of input matrix M */
  DEM           = (double *)mxGetPr(prhs[0]); 
  M             = mxGetM(prhs[0]);
  N             = mxGetN(prhs[0]);
  Pud_pts       = (double *)mxGetPr(prhs[1]); 
  Mp            = mxGetM(prhs[1]);

  Thres_pts     = (double *)mxGetPr(prhs[2]);   
  Mt            = mxGetM(prhs[2]);

  ind_pudcum    = (double *)mxGetPr(prhs[3]); 
  num_pud       = mxGetNumberOfElements(prhs[3])-1;

  ind_threscum  = (double *)mxGetPr(prhs[4]); 
  num_thres     = mxGetNumberOfElements(prhs[4])-1;

  Center_real   = (double *)mxGetPr(prhs[5]);
  dx            = mxGetScalar(prhs[6]);
  dy            = mxGetScalar(prhs[7]);

  volume        = (double*)malloc(num_pud*sizeof(double));  /* Initialize volume */
  depth         = (double*)malloc(num_pud*sizeof(double));  /* Initialize depth */

  for (i=0;i<num_pud;i++){
    volume[i] = 0;
    k = (int)ind_threscum[i];
    depth[i] = DEM[M*(int)Thres_pts[0*Mt+k]+(int)Thres_pts[1*Mt+k]] - Center_real[i];

    for (j=(int)ind_pudcum[i];j<(int)ind_pudcum[i+1];j++){
      volume[i] += (DEM[M*(int)Thres_pts[0*Mt+k]+(int)Thres_pts[1*Mt+k]]-DEM[M*(int)Pud_pts[0*Mp+j]+(int)Pud_pts[1*Mp+j]])*dx*dy;
    }
  }

  plhs[0] = mxCreateNumericMatrix(num_pud,1, mxDOUBLE_CLASS, mxREAL);
  memmove(mxGetData(plhs[0]),volume,num_pud*sizeof(double));

  plhs[1] = mxCreateNumericMatrix(num_pud,1, mxDOUBLE_CLASS, mxREAL);
  memmove(mxGetData(plhs[1]),depth,num_pud*sizeof(double));

} /* end mexFunction() */
