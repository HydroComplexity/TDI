/*=======================================================================
* mex_DEM_correct.c - C MEX file used to fill all gaps and missing data 
      of multi centers to their lowest neighbor
*------------------------------------------------------------------------
*     usage: mex_DEM_correct_mul(DEM, Pud_pts, ind_pudcum)

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
  this copyright notice is not removed. 
*========================================================================
*/

#include "mex.h"
#include <string.h>
#include "matrix.h"
#include <math.h>

/* ------ FUNCTION PROPOTYPES ---------------------------------------------*/
void getmin(double MAT[], int M, int N, int arr_in[], int size_in, double *min_out);
void getdiff(int arr1[], int size1, int arr2[], int size2, int arr_out[], int *num_out);
void getneigh(int arr_size, int arr[], int arr_out[]);


/* ------ BEGIN mexFunction() ---------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *DEM, *Psi, *Cen_pts, *ind_pts, *ind_cum, nodata, minval;
  int *Af, *Af2, *Df, *Diff;
  int i, j, p, M, N, Mp, num_pud, num_diff, mink;
  mxArray *output_1, *output_2[2], *input_1, *input_2[2];

  /* Get value and dimensions of input matrix M */
  DEM     = (double *)mxGetPr(prhs[0]); 
  M       = mxGetM(prhs[0]);
  N       = mxGetN(prhs[0]);

  Cen_pts = (double *)mxGetPr(prhs[1]);
  Mp      = mxGetN(prhs[1]);
  ind_pts = (double *)mxGetPr(prhs[2]);
  ind_cum = (double *)mxGetPr(prhs[3]);
  num_pud = mxGetNumberOfElements(prhs[3]);
  nodata  = mxGetScalar(prhs[4]);
  mink    = mxGetScalar(prhs[4]);

  Af      = (int*)malloc(2*sizeof(int));
  Df      = (int*)malloc(sizeof(int));
  Diff    = (int*)malloc(sizeof(int));
  Psi     = (double*)malloc(M*N*sizeof(double));
  memcpy(Psi, DEM, M*N*sizeof(double));	

  for (i=1;i<num_pud;i++){
    if (Psi[M*(int)Cen_pts[2*(int)ind_cum[i-1]]+(int)Cen_pts[2*(int)ind_cum[i-1]+1]] == nodata){
      p   = (int)ind_pts[i];
      Af  = realloc(Af,2*p*sizeof(int));
      for (j=(int)ind_cum[i-1];j<(int)ind_cum[i];j++){
        Af[2*(j-(int)ind_cum[i-1])]   = (int)Cen_pts[2*j];
        Af[2*(j-(int)ind_cum[i-1])+1] = (int)Cen_pts[2*j+1];
      }
      if (p<mink){
        Df = realloc(Df,16*p*sizeof(int));   
        getneigh(p, Af, Df);

        Diff = realloc(Diff,16*p*sizeof(int));
        getdiff(Df, 8*p, Af, p, Diff, &num_diff);
      }
      else{
        input_1 = mxCreateNumericMatrix(p,2,mxINT32_CLASS,mxREAL);
        Af2     = mxGetData(input_1);	
        memcpy(Af2,Af,p*2*sizeof(int));
        mexCallMATLAB(1,&output_1,1,&input_1,"minksumMEX");
        input_2[0] = output_1;	
        input_2[1] = input_1;

        mexCallMATLAB(2,output_2,2,input_2,"matdiff_mex");
        num_diff  = mxGetScalar(output_2[1]);
        Diff      = realloc(Diff,2*num_diff*sizeof(int));   
        memcpy(Diff,mxGetData(output_2[0]), 2*num_diff*sizeof(int));
      }
      getmin(DEM, M, N, Diff, 8, &minval);

      for (j=(int)ind_cum[i-1];j<(int)ind_cum[i];j++){
        Psi[M*(int)Cen_pts[2*j]+(int)Cen_pts[2*j+1]] = minval-0.01;
      }
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
/* ------ Get Difference Cells---------------------------------------------
  Get different cells of array 1 and array 2
  Getdiff will return the cells surround the old depression cells
---------------------------------------------------------------------------*/
void getdiff(int arr1[], int size1, int arr2[], int size2, int arr_out[], int *num_out)
{
  int i, j, p, q, K, flag;
  int *arr_temp;

  K = 2;
  p = 0;
  q = 0;

  arr_temp  = (int*)malloc(2*size1*sizeof(int));

  for (i=0;i<size1;i++){
    flag = 1;
    for (j=0;j<size2;j++){
      if (arr1[K*i] == arr2[K*j] && arr1[K*i+1] == arr2[K*j+1]) {
        flag = 0;
      }
    }
    if (flag ==1){
      arr_temp[K*p]   = arr1[K*i];
      arr_temp[K*p+1] = arr1[K*i+1];
      p++;
    }
  }
  arr_temp 	= realloc(arr_temp,2*p*sizeof(int));

  for (i=0;i<p-1;i++){
    for (j=i+1;j<p;j++){
      if (arr_temp[K*i]!= -999 && arr_temp[K*i+1] != -999 && arr_temp[K*i]==arr_temp[K*j] && arr_temp[K*i+1]==arr_temp[K*j+1]){
        arr_temp[K*j]   = -999;
        arr_temp[K*j+1] = -999;
      }
    }
  }
  for (i=0;i<p;i++){
    if (arr_temp[K*i] != -999 && arr_temp[K*i+1] != -999){
      arr_out[K*q]    = arr_temp[K*i];
      arr_out[K*q+1]  = arr_temp[K*i+1];
      q++;
    }
  }
  arr_out 	= realloc(arr_out,2*q*sizeof(int));
  *num_out	= q;
}


/* ------ Get Neighbor Cells---------------------------------------------
  Get neighbor cells around point(s)
  D8 neighbor of all cellls
-------------------------------------------------------------------------*/
void getneigh(int arr_size, int arr[], int arr_out[])
{
  int i, j, p, q, K, flag;

  arr_out = realloc(arr_out,16*arr_size*sizeof(int));
  K       = 2;
  p       = 0;

  for (q=0;q<arr_size;q++){
    arr_out[K*p+0]      = arr[K*q+0] - 1;		
    arr_out[K*p+1]      = arr[K*q+1] - 1;
    arr_out[K*(p+1)+0]  = arr[K*q+0] - 1;		
    arr_out[K*(p+1)+1]  = arr[K*q+1];

    arr_out[K*(p+2)+0]  = arr[K*q+0] - 1;		
    arr_out[K*(p+2)+1]  = arr[K*q+1] + 1;
    arr_out[K*(p+3)+0]  = arr[K*q+0];   
    arr_out[K*(p+3)+1]  = arr[K*q+1] - 1;

    arr_out[K*(p+4)+0]  = arr[K*q+0];		
    arr_out[K*(p+4)+1]  = arr[K*q+1] + 1;
    arr_out[K*(p+5)+0]  = arr[K*q+0] + 1;		
    arr_out[K*(p+5)+1]  = arr[K*q+1] - 1;

    arr_out[K*(p+6)+0]  = arr[K*q+0] + 1;		
    arr_out[K*(p+6)+1]  = arr[K*q+1];
    arr_out[K*(p+7)+0]  = arr[K*q+0] + 1;		
    arr_out[K*(p+7)+1]  = arr[K*q+1] + 1;
    p += 8;
  }
}

/* ------ END OF ALL FUNCTION ------------------------------------------*/
