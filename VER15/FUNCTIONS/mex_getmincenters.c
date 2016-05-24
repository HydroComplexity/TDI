/*=======================================================================
* mex_getmincenter.c - C MEX file used to identify the true Centers
    of depression cells in high levels due to DEM filling processes

  In levels 2,3,4..n, the cells are filled with threshold.
  We need to identified the lowest cell of flat in such level for
  accurate depth and volume.		
*------------------------------------------------------------------------
*     usage: [Center] = mex_getmincenter(DEM,Cen_multi,);

*------------------------------------------------------------------------
*     Input arguments:
*       DEM:          a K x J matrix of elevation
*       Pud_pts:      List of depression cells
*       Thres_pts:    List of depression cells
*       ind_pudcum:   Cumulative index of ind_pud
*       ind_threscum: Cumulative index of ind_thres
*
*     Return arguments:
*       Volume:       list of elevation of real centers 
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
  double *DEM, *Cen_multi, *ind_mulcum, *Cen_flats;
  int i, j, M, N, k, num_cenmul, Mc;
  double min_val;

  /* Get value and dimensions of input matrix M */
  DEM         = (double *)mxGetPr(prhs[0]);
  M           = mxGetM(prhs[0]);
  N           = mxGetN(prhs[0]);
  Cen_multi   = (double *)mxGetPr(prhs[1]);
  Mc          = mxGetM(prhs[1]);
  ind_mulcum  = (double *)mxGetPr(prhs[2]);
  num_cenmul  = mxGetNumberOfElements(prhs[2])-1;
  Cen_flats   = (double*)malloc(num_cenmul*sizeof(double)); /* Initialize Center flats */

  for (i=0;i<num_cenmul;i++){
    min_val = DEM[M*(int)Cen_multi[0*Mc+(int)ind_mulcum[i]]+(int)Cen_multi[1*Mc+(int)ind_mulcum[i]]];
    for (j=(int)ind_mulcum[i];j<(int)ind_mulcum[i+1];j++){
      if (DEM[M*(int)Cen_multi[0*Mc+j]+(int)Cen_multi[1*Mc+j]] < min_val){
        min_val = DEM[M*(int)Cen_multi[0*Mc+j]+(int)Cen_multi[1*Mc+j]];
      }
    }
    Cen_flats[i] = min_val;
  }
  plhs[0] = mxCreateNumericMatrix(num_cenmul,1, mxDOUBLE_CLASS, mxREAL);
  memmove(mxGetData(plhs[0]),Cen_flats,num_cenmul*sizeof(double));	
}
