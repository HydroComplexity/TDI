/*=======================================================================
* mex_fill_threshold.c - C MEX file used to fill all depression cells to 
      their threshold level value
*------------------------------------------------------------------------
* usage: mex_fill_threshold(Psi,Pud_pts,Thres_pts,ind_pudcum,ind_threscum,ind_one);

*------------------------------------------------------------------------
*     Input arguments:
*       DEM:          a M x N matrix of elevation
*       Pud_pts:      List of depression cells
*       Thres_pts:    List of depression cells
*       ind_pudcum:   Cumulative index of ind_pud
*       ind_threscum: Cumulative index of ind_thres
*       ind_one:      index of ONE matrix from Matlab
*
*     Return arguments:
*       No return argument, the function directly changes the pointer
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

/* ------ FUNCTION PROPOTYPES ---------------------------------------------*/
void getbounds(double arr_in[], int size_in, double num[], int num_one, int *lowb, int *upb);


/* ------ BEGIN mexFunction() ---------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *DEM, *Psi, *Pud_pts, *Thres_pts, *ind_pudcum, *ind_threscum, *ind_one;
  int i, j, M, N, Mp, Mt, num_thres, num_one, num_pud, *lowb, *upb;

  /* Get value and dimensions of input matrix M */
  DEM           = (double *)mxGetPr(prhs[0]); 
  M             = mxGetM(prhs[0]);
  N             = mxGetN(prhs[0]);

  Pud_pts       = (double *)mxGetPr(prhs[1]);
  Mp            = mxGetM(prhs[1]);

  Thres_pts     = (double *)mxGetPr(prhs[2]);
  Mt            = mxGetM(prhs[2]);

  ind_pudcum    = (double *)mxGetPr(prhs[3]);

  ind_threscum  = (double *)mxGetPr(prhs[4]);
  num_thres     = mxGetNumberOfElements(prhs[4]);

  ind_one       = (double *)mxGetPr(prhs[5]);
  num_one       = mxGetNumberOfElements(prhs[5]);

  lowb          = (int*)malloc(num_one*sizeof(int));
  upb           = (int*)malloc(num_one*sizeof(int));

  getbounds(ind_threscum, num_thres, ind_one, num_one, lowb, upb);

  for (i=0;i<num_one;i++){
    for (j=(int)ind_pudcum[lowb[i]]; j<(int)ind_pudcum[upb[i]];j++){
      DEM[M*(int)Pud_pts[0*Mp+j]+(int)Pud_pts[1*Mp+j]] = DEM[M*(int)Thres_pts[0*Mt+(int)ind_one[i]]+(int)Thres_pts[1*Mt+(int)ind_one[i]]];
    }
  }
} /* end mexFunction() */

/* ------ Get Boundary Cells-------------------------------------------------
---------------------------------------------------------------------------*/
void getbounds(double arr_in[], int size_in, double num[], int num_one, int *lowb, int *upb)
{
  int i,j, k, jf;
  j = 0;
  while (j<num_one){
    for (i=0;i<size_in;i++){
      if ((int) arr_in[i] > (int)num[j] && j<num_one){
        lowb[j] = i-1;
        upb[j]  = i;
        j++;
        jf = j;

        for (k=jf;k<num_one;k++){
          if ((int)arr_in[i] > (int)num[k]){
            lowb[j] = i-1;
            upb[j]  = i;
            j++;
          }
          else{
            break;
          }
        }
      }
    }
  }
}
