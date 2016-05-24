/*=========================================================================
* mex_mul_centhres.c - C MEX file used to get the locations of Multi-center
*     cells from input image M.
*--------------------------------------------------------------------------
*     usage: [A_multi] = 
        mex_mul_centhres(DEM,S,linC,len_grenum,ind_hist,ind_mul,Bs)
*--------------------------------------------------------------------------
*     Input arguments:
*       DEM:        a M x N matrix of elevations
*       S:          Sparse matrix of connectivity image B from Matlab
*       linC:       unique value of S
*       len_grenum: length of ind_hist (Matlab) greater than 1
*       ind_hist:   index histogram from mex_regHistogram
*       ind_mul:    index multi center		
*
*     Return arguments:
*       A_multi:    a array (n x 2) of zero-value (multi-center) cells
*--------------------------------------------------------------------------
* This is a MEX-file for MATLAB.
* Created by: Phong Le <phongle1@illinois.edu>
* This code may be freely redistributed in its entirety provided that 
  this copyright notice is	not removed.
*==========================================================================
*/

#include "mex.h"
#include <string.h>
#include "matrix.h"
#include <math.h>

/* ------ FUNCTION PROPOTYPES ---------------------------------------------*/

void getdiff(int arr1[], int size1, int arr2[], int size2, int arr_out[], int *num_out);
void getneigh(int arr_size, int arr[], int arr_out[]);


/* ------ BEGIN mexFunction() ---------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize M, N, nz, nrow;
  double *pr, *DEM, *linC, *ind_hist, *ind_mul, *Bs, elev;
  mwIndex *ir, *jc;
  int i, j, k, h, p, t,num_neigh,flag_neigh,flag_cen,len_grenum,Amulsize,mink;
  int *Af, *Af2, *Afa, *Df, *neighbor, *A_multi;
  mxArray *output_1, *output_2[2], *input_1, *input_2[2];

  if (! mxIsSparse (prhs[1]))
    mexErrMsgTxt ("expects sparse matrix");

  DEM         = (double *)mxGetPr(prhs[0]); 														
  M           = mxGetM (prhs [0]);
  N           = mxGetN (prhs [0]);

  pr          = mxGetPr (prhs[1]);
  ir          = mxGetIr (prhs[1]);
  jc          = mxGetJc (prhs[1]);

  linC        = (double*)mxGetPr(prhs[2]);

  len_grenum  = mxGetScalar(prhs[3]);
  ind_hist    = (double*)mxGetPr(prhs[4]);
  ind_mul     = (double*)mxGetPr(prhs[5]);
  Bs          = (double*)mxGetPr(prhs[6]);
  mink        = mxGetScalar(prhs[7]);

  neighbor    = (int*)mxMalloc(sizeof(int));
  Af          = (int*)malloc(2*sizeof(int));
  Af2         = (int*)malloc(2*sizeof(int));
  Df          = (int*)mxMalloc(sizeof(int));
  Afa         = (int*)malloc(sizeof(int));
  A_multi     = (int*)malloc(2*sizeof(int));
  A_multi[0]  = -1000;
  A_multi[1]  = -1000;
  Amulsize    = 1;

  for (k=1; k<len_grenum;k++){
    h = 0;
    p = 0;
    for (j=0;j<N;j++) {
      nrow = jc[j+1] - jc[j];
      for (i=0;i<nrow;i++){	
        if (pr[h] == linC[(int)ind_hist[k]-1]){
          Af        = realloc(Af,2*(p+1)*sizeof(int));
          Af[2*p]   = j;
          Af[2*p+1] = ir[h];
          p++;
        }
        h++;
      }
    }
    /* If number of point is >= mink, use MATLAB minkowski function	*/		
    if (p<mink){
      Df = mxRealloc(Df,16*p*sizeof(int));
      getneigh(p,Af,Df);

      neighbor = mxRealloc(neighbor,16*p*sizeof(int));
      getdiff(Df,8*p,Af,p,neighbor,&num_neigh);
    }
    else{
      input_1     = mxCreateNumericMatrix(p, 2, mxINT32_CLASS, mxREAL);
      Af2         = mxGetData(input_1);	
      memcpy(Af2,Af,p*2*sizeof(int));
      mexCallMATLAB(1,&output_1,1,&input_1,"minksumMEX");
      input_2[0]  = output_1;	
      input_2[1]  = input_1;

      mexCallMATLAB(2,output_2,2,input_2,"matdiff_mex");
      num_neigh = mxGetScalar(output_2[1]);
      neighbor  = mxRealloc(neighbor,2*num_neigh*sizeof(int));   
      memcpy(neighbor,mxGetData(output_2[0]), 2*num_neigh*sizeof(int));
    }
    flag_neigh = 0;
    for (j=0;j<num_neigh;j++){
      if (neighbor[2*j]<0 || neighbor[2*j]>N-1 || neighbor[2*j+1]<0 || neighbor[2*j+1]>M-1){
        flag_neigh = 1;
        break;
      }
    }
    if (flag_neigh == 0){
      elev = DEM[M*Af[0]+Af[1]];
      t = 0;
      for (j=0;j<num_neigh;j++){
        if (DEM[M*neighbor[2*j]+neighbor[2*j+1]]==elev){
          Afa         = realloc(Afa,2*(t+1)*sizeof(int));
          Afa[2*t]    = neighbor[2*j];
          Afa[2*t+1]  = neighbor[2*j+1];
          t++; 
        }
      }
      Af = realloc(Af,2*(p+t)*sizeof(int)); 
      memcpy(&Af[2*p],Afa,2*t*sizeof(int)); 
      
      if (p+t<mink){
        Df = mxRealloc(Df,16*(p+t)*sizeof(int));   
        getneigh(p+t,Af,Df);
        neighbor = mxRealloc(neighbor,16*(p+t)*sizeof(int));
        getdiff(Df,8*(p+t),Af,p,neighbor,&num_neigh);
      }
      else{
        input_1     = mxCreateNumericMatrix(p+t, 2, mxINT32_CLASS, mxREAL);
        Af2         = mxGetData(input_1);
        memcpy(Af2,Af,(p+t)*2*sizeof(int));
        mexCallMATLAB(1,&output_1,1,&input_1,"minksumMEX");
        input_2[0]  = output_1;	
        input_2[1]  = input_1;

        mexCallMATLAB(2,output_2,2,input_2,"matdiff_mex");
        num_neigh   = mxGetScalar(output_2[1]);
        neighbor    = mxRealloc(neighbor,2*num_neigh*sizeof(int));   
        memcpy(neighbor,mxGetData(output_2[0]), 2*num_neigh*sizeof(int));
      }
    }

    for (j=0;j<p;j++)
      Bs[M*Af[2*j]+Af[2*j+1]] = 0;

    flag_cen = 0;
    for (j=0;j<num_neigh;j++){
      if (neighbor[2*j] >= 0 && neighbor[2*j] < N && neighbor[2*j+1] >=0 && neighbor[2*j+1] < M){
        if (DEM[M*neighbor[2*j]+neighbor[2*j+1]] <= DEM[M*Af[0]+Af[1]])
          flag_cen = 1;
      }	
      else{
        flag_cen = 1;
      }
    }
    if (flag_cen==0){
      A_multi = realloc(A_multi,2*(Amulsize+p+t)*sizeof(int));  
      memcpy(&A_multi[2*Amulsize],Af,2*(p+t)*sizeof(int));
      Amulsize += p+t;
      ind_mul[k] = (double)(p+t);
    }
  }
  plhs[0] = mxCreateNumericMatrix(2,Amulsize, mxINT32_CLASS, mxREAL);
  memmove(mxGetData(plhs[0]),A_multi,2*Amulsize*sizeof(int));
} /* END mexFunction() */


/* ------ Get Difference Cells---------------------------------------------
  Get different cells of array 1 and array 2                                
---------------------------------------------------------------------------*/
void getdiff(int arr1[], int size1, int arr2[], int size2, int arr_out[], int *num_out)
{
  int i, j, p, q, K, flag;
  int *arr_temp;

  K = 2;
  p = 0;
  q = 0;

  arr_temp	= (int*)malloc(2*size1*sizeof(int));
  for (i=0;i<size1;i++){
    flag = 1;
    for (j=0;j<size2;j++){
      if (arr1[K*i] == arr2[K*j] && arr1[K*i+1] == arr2[K*j+1]) {
        flag = 0;
        break;
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
      if (arr_temp[K*i] != -999 && arr_temp[K*i+1] != -999 && arr_temp[K*i]==arr_temp[K*j] && arr_temp[K*i+1]==arr_temp[K*j+1]){
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
  arr_out   = mxRealloc(arr_out,2*q*sizeof(int));
  *num_out  = q;
}


/* ------ Get Neighbor Cells---------------------------------------------
  Get neighbor cells around point(s)
-------------------------------------------------------------------------*/
void getneigh(int arr_size, int arr[], int arr_out[])
{
  int i, j, p, q, K, flag;

  K = 2;
  p = 0;
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
