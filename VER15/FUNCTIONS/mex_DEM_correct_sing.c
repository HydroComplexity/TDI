/*=======================================================================
* mex_DEM_correct_sing.c - C MEX file used to fill all gaps and missing data 
      of single center to their lowest neighbor
*------------------------------------------------------------------------
*     usage: mex_DEM_correct_sing(DEM, Pud_pts, ind_pudcum)

*------------------------------------------------------------------------
*     Input arguments:
*       DEM:     		a M x N matrix of elevation
*       Pud_pts:			List of depression cells
*       ind_pudcum:	Cumulative index of ind_pud
*
*     Return arguments:
*       DEM_pud     	a new M x N matrix of elevation
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
void getmin(double MAT[], int M, int N, int arr_in[], int size_in, double *min_out);
void getneigh(int arr[], int arr_out[]);


/* ------ BEGIN mexFunction() ---------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *DEM, *Psi, *Cen_pts, nodata, minval;
  int *Af, *Df;
  int i, j, p, M, N, num_pud;

  /* Get value and dimensions of input matrix M */
  DEM     = (double *)mxGetPr(prhs[0]); 														
  M       = mxGetM(prhs[0]);
  N       = mxGetN(prhs[0]);	

  Cen_pts = (double *)mxGetPr(prhs[1]); 														
  num_pud = mxGetN(prhs[1]);
  nodata  = mxGetScalar(prhs[2]);	

  Af      = (int*)malloc(2*sizeof(int));
  Df      = (int*)malloc(16*sizeof(int));
  Psi     = (double*)malloc(M*N*sizeof(double));
  memcpy(Psi, DEM, M*N*sizeof(double));	

  for (i=0;i<num_pud;i++){
    if (Psi[M*(int)Cen_pts[2*i]+(int)Cen_pts[2*i+1]] == nodata){
      Af[0] = (int)Cen_pts[2*i];
      Af[1] = (int)Cen_pts[2*i+1];
      getneigh(Af, Df);

      getmin(DEM, M, N, Df, 8, &minval);
      Psi[M*(int)Cen_pts[2*i]+(int)Cen_pts[2*i+1]] = minval-0.01;
    } 
  }

  plhs[0] = mxCreateNumericMatrix(M,N, mxDOUBLE_CLASS, mxREAL);
  memmove(mxGetData(plhs[0]),Psi,M*N*sizeof(double));

} /* end mexFunction() */


/* ------ Get Minimum Value ---------------------------------------------
  Get Min values of array_in, and return to one value min_out.
-------------------------------------------------------------------------*/
void getmin(double MAT[], int M, int N, int arr_in[], int size_in, double *min_out)
{
  int i,  p, K;
  double min_val;
  K = 2;
  p = 0;

  min_val = MAT[M*arr_in[0]+arr_in[1]];
  for (i=1;i<size_in;i++){
    if (MAT[M*arr_in[K*i]+arr_in[K*i+1]]<min_val){
      min_val = MAT[M*arr_in[K*i]+arr_in[K*i+1]];
    }
  }
  *min_out = min_val;
}

/* ------ Get Neighbor Cells---------------------------------------------
  Get neighbor cells around point(s)
  D8 neighbor of all cellls
-------------------------------------------------------------------------*/
void getneigh(int arr[], int arr_out[])
{
  int p;
  p = 0;

  arr_out[2*p+0]      = arr[0] - 1;
  arr_out[2*p+1]      = arr[1] - 1;
  arr_out[2*(p+1)+0]  = arr[0] - 1;
  arr_out[2*(p+1)+1]  = arr[1];

  arr_out[2*(p+2)+0]  = arr[0] - 1;
  arr_out[2*(p+2)+1]  = arr[1] + 1;
  arr_out[2*(p+3)+0]  = arr[0];		
  arr_out[2*(p+3)+1]  = arr[1] - 1;

  arr_out[2*(p+4)+0]  = arr[0];		
  arr_out[2*(p+4)+1]  = arr[1] + 1;
  arr_out[2*(p+5)+0]  = arr[0] + 1;
  arr_out[2*(p+5)+1]  = arr[1] - 1;

  arr_out[2*(p+6)+0]  = arr[0] + 1;
  arr_out[2*(p+6)+1]  = arr[1];
  arr_out[2*(p+7)+0]  = arr[0] + 1;
  arr_out[2*(p+7)+1]  = arr[1] + 1;
}

/* ------ END OF ALL FUNCTION ------------------------------------------*/
