/*=========================================================================
* mex_other_run.c - C MEX file used to get the locations of depression
*	      cells from input image M and cen_multi list.
*--------------------------------------------------------------------------
*     usage: [Pud_pts, Thres_pts,ind_pud, ind_thres] = 
        mex_others_run(MAT,Cen_pts,num_pts,ind_pts,ind_cum,Puddle,Threshold);
*--------------------------------------------------------------------------
*     Input arguments:
*       MAT:        a M x N matrix of elevations
*       S:          Sparse matrix of connectivity image B from Matlab
*       Cen_pts:    List of Center cell points
*       num_pts:    Number of points
*       nd_pts:     Index of Cen_pts
*       ind_cum:    Cumulative index of ind_pts
*       Puddle:     a M x N matrix of Puddle cells
*       Threshold:  a M x N matrix of Threshold cells
*
*     Return arguments:
*       Pud_pts:    List of Puddle cell points
*       Thres_pts:  List of Threshold cell points
*       ind_pts:    Index of Puddle cells
*       ind_pts:    Index of Threshold cells
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


/* ------ FUNCTION PROTOTYPES ---------------------------------------------*/
void getmin(double MAT[], int M, int N, int arr_in[], int size_in, int arr_out[], int *num_out);
void getdiff(int arr1[], int size1, int arr2[], int size2, int arr_out[], int *num_out);
void getneigh(int arr_size, int arr[], int arr_out[]);

/* ------ BEGIN mexFunction() ---------------------------------------------*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *Psi, *DEM, *Cen_pts, *ind_pts, *ind_cum, *Thres_in, *Puddle_in, *Puddle, *Threshold;	
  int *Apud, *Apud2, *Dpud, *Tpud, *Diffpud, *Diff_min, *min_list,*Dlist_min, *Sub_min;
  int *Pud_pts, *Thres_pts, *ind_pud, *ind_thres;
  int M, N, loop_sear, thressize, pudsize, num_pts, Cen_size, num_Dpud, num_diff, num_min, num_submin, count;
  int sum_thressize, sum_pudsize, sum_thressize_old, sum_pudsize_old;  
  int i,j,k,l,mink,max_size;
  mxArray *output_1, *output_2[2], *input_1, *input_2[2];

  /* Get value and dimensions of input matrix M */
  DEM         = (double *)mxGetPr(prhs[0]);
  M           = mxGetM(prhs[0]);
  N           = mxGetN(prhs[0]);
  Cen_pts     = (double *)mxGetPr(prhs[1]);
  num_pts     = mxGetScalar(prhs[2]);/* num */
  ind_pts     = mxGetPr(prhs[3]);/* length num+1  */
  ind_cum     = mxGetPr(prhs[4]); /*	length num+1	*/
  Puddle      = (double *)mxGetPr(prhs[5]); 
  Threshold   = (double *)mxGetPr(prhs[6]); 
  mink        = mxGetScalar(prhs[7]); 
  max_size    = mxGetScalar(prhs[8]); 

  Cen_size    = mxGetN(prhs[1]); /*Total number of center points*/
  Apud        = (int*)malloc(sizeof(int));
  Dpud        = (int*)malloc(sizeof(int));
  Pud_pts     = (int*)malloc(2*sizeof(int));
  Thres_pts   = (int*)malloc(2*sizeof(int));
  Tpud        = (int*)malloc(2*sizeof(int));

  Diffpud     = (int*)mxMalloc(sizeof(int));
  Diff_min    = (int*)mxMalloc(sizeof(int));
  Sub_min     = (int*)mxMalloc(sizeof(int));

  ind_pud     = (int*)malloc((num_pts+1)*sizeof(int));
  ind_thres   = (int*)malloc((num_pts+1)*sizeof(int));

  min_list    = (int*)malloc(2*sizeof(int));
  Dlist_min   = (int*)malloc(16*sizeof(int));
  Psi         = (double*)malloc(M*N*sizeof(double));
  memcpy(Psi, DEM, M*N*sizeof(double));

  ind_pud[0]        = 1;
  ind_thres[0]      = 1;
  sum_pudsize_old   = 0;
  sum_thressize_old = 0;
  sum_pudsize       = 0;
  sum_thressize     = 0;

  for (k=0;k<num_pts;k++){
    pudsize = (int)ind_pts[k+1];
    Apud 		= realloc(Apud,2*pudsize*sizeof(int));  /* re-allocate storage for the position center cells  */

  	for (j=0;j<pudsize;j++){
      for (i=0;i<2;i++){
        Apud[2*j+i] = (int) Cen_pts[2*(int)ind_cum[k]+2*j+i];
      }
      Psi[M*Apud[2*j] + Apud[2*j+1]]    = 99999;
      Puddle[M*Apud[2*j] + Apud[2*j+1]] = 1.0;
    }

    if (pudsize < max_size){
      loop_sear = 0;
      count     = 0;
      while (loop_sear<1 && count < 2*max_size){
        if (pudsize < mink){
          Dpud = realloc(Dpud,16*pudsize*sizeof(int));   
          getneigh(pudsize, Apud, Dpud);

          Diffpud = mxRealloc(Diffpud,16*pudsize*sizeof(int));
          getdiff(Dpud, 8*pudsize, Apud, pudsize, Diffpud, &num_diff);
        }
        else{
          input_1 = mxCreateNumericMatrix(pudsize, 2, mxINT32_CLASS, mxREAL);
          Apud2   = mxGetData(input_1);
          memcpy(Apud2,Apud,pudsize*2*sizeof(int));
          mexCallMATLAB(1,&output_1,1,&input_1,"minksumMEX");
          input_2[0] = output_1;
          input_2[1] = input_1;

          mexCallMATLAB(2,output_2,2,input_2,"matdiff_mex");
          num_diff  = mxGetScalar(output_2[1]);
          Diffpud   = mxRealloc(Diffpud,2*num_diff*sizeof(int));   
          memcpy(Diffpud,mxGetData(output_2[0]), 2*num_diff*sizeof(int));
        }
        Diff_min = mxRealloc(Diff_min,2*num_diff*sizeof(int));
        getmin(Psi, M, N, Diffpud, num_diff, Diff_min, &num_min);

        thressize = 0;
        for (i=0;i<num_min;i++){
          min_list[0] = Diff_min[2*i];
          min_list[1] = Diff_min[2*i+1];
          getneigh(1, min_list, Dlist_min);

          Sub_min = mxRealloc(Sub_min,16*sizeof(int));
          getmin(Psi, M, N, Dlist_min, 8, Sub_min, &num_submin);

          if (Psi[M*min_list[0]+min_list[1]] < Psi[M*Sub_min[0]+Sub_min[1]] && min_list[0]>0 && min_list[0] < N-1 && min_list[1]>0 && min_list[1]<M-1){
            Puddle[M*min_list[0] + min_list[1]] = 1.0;
            Psi[M*min_list[0] + min_list[1]]    = 99999;
            pudsize++;
            Apud = realloc(Apud,2*pudsize*sizeof(int));
            Apud[2*(pudsize-1)]   = min_list[0];
            Apud[2*(pudsize-1)+1] = min_list[1];
          }
          else{
            thressize++;
            Threshold[M*min_list[0]+min_list[1]] = Threshold[M*min_list[0]+min_list[1]] + 1;

            Tpud = realloc(Tpud,2*thressize*sizeof(int));
            Tpud[2*(thressize-1)]   = min_list[0];
            Tpud[2*(thressize-1)+1] = min_list[1];
            loop_sear	= 1;
          }
        }
        count++;
      }
    }
    else  {
      input_1 = mxCreateNumericMatrix(pudsize, 2, mxINT32_CLASS, mxREAL);
      Apud2   = mxGetData(input_1);	
      memcpy(Apud2,Apud,pudsize*2*sizeof(int));
      mexCallMATLAB(1,&output_1,1,&input_1,"minksumMEX");
      input_2[0] = output_1;
      input_2[1] = input_1;

      mexCallMATLAB(2,output_2,2,input_2,"matdiff_mex");
      num_diff  = mxGetScalar(output_2[1]);
      Diffpud   = mxRealloc(Diffpud,2*num_diff*sizeof(int));   
      memcpy(Diffpud,mxGetData(output_2[0]), 2*num_diff*sizeof(int));

      Diff_min = mxRealloc(Diff_min,2*num_diff*sizeof(int));
      getmin(Psi, M, N, Diffpud, num_diff, Diff_min, &num_min);

      thressize = num_min;
      Tpud = realloc(Tpud,2*num_min*sizeof(int));
      memcpy(Tpud,Diff_min,2*num_min*sizeof(int));
    }	
    ind_pud[k+1]    = pudsize;
    ind_thres[k+1]  = thressize;
    sum_pudsize     += pudsize;
    sum_thressize   += thressize;

    Pud_pts     = realloc(Pud_pts, 2*sum_pudsize*sizeof(int));  
    Thres_pts   = realloc(Thres_pts, 2*sum_thressize*sizeof(int));

    memcpy(&Pud_pts[2*sum_pudsize_old], Apud, 2*pudsize*sizeof(int));
    memcpy(&Thres_pts[2*sum_thressize_old], Tpud, 2*thressize*sizeof(int));

    sum_pudsize_old   = sum_pudsize;
    sum_thressize_old = sum_thressize;
    /* Copy back the elevation for searching      */
    for (i=0;i<pudsize;i++){
      Psi[M*Apud[2*i]+Apud[2*i+1]] = DEM[M*Apud[2*i]+Apud[2*i+1]];
    }
  }	
  /* Return array & matrix . . . . */
  plhs[0] = mxCreateNumericMatrix(sum_pudsize,2, mxINT32_CLASS, mxREAL);
  memmove(mxGetData(plhs[0]),Pud_pts,2*sum_pudsize*sizeof(int));

  plhs[1] = mxCreateNumericMatrix(sum_thressize,2, mxINT32_CLASS, mxREAL);
  memmove(mxGetData(plhs[1]),Thres_pts,2*sum_thressize*sizeof(int));

  plhs[2] = mxCreateNumericMatrix(num_pts+1,1, mxINT32_CLASS, mxREAL);
  memmove(mxGetData(plhs[2]),ind_pud,(num_pts+1)*sizeof(int));

  plhs[3] = mxCreateNumericMatrix(num_pts+1,1, mxINT32_CLASS, mxREAL);
  memmove(mxGetData(plhs[3]),ind_thres,(num_pts+1)*sizeof(int));

} /* end mexFunction() */



/* ------ Get Minimum Value ---------------------------------------------
  Get Min values of array_in, and return to arr_out.
  If the set has > 1 min, the function returns a vector
-------------------------------------------------------------------------*/
void getmin(double MAT[], int M, int N, int arr_in[], int size_in, int arr_out[], int *num_out)
{
  int i, j, p, K;
  double min_val;

  K = 2;
  p = 0;
  min_val = 99999;

  for (i=0;i<size_in;i++){
    if (arr_in[K*i]>=0 && arr_in[K*i]<N && arr_in[K*i+1]>=0 && arr_in[K*i+1]<M && MAT[M*arr_in[K*i]+arr_in[K*i+1]]<min_val) {
      min_val = MAT[M*arr_in[K*i]+arr_in[K*i+1]];
    }
  }

  for (i=0;i<size_in;i++){
    if (arr_in[K*i]>=0 && arr_in[K*i]<N && arr_in[K*i+1]>=0 && arr_in[K*i+1]<M && MAT[M*arr_in[K*i]+arr_in[K*i+1]]==min_val) {
      arr_out[K*p] 		= arr_in[K*i];
      arr_out[K*p+1] 	= arr_in[K*i+1];
      p++;
    }
  }

  if (p==0){
    printf("Warning: p = %d, min_val = %f6.3 and MAT = %f6.3\n",p, min_val,MAT[M*arr_in[K*0]+arr_in[K*0+1]]);
    printf("Check the min_val default parameters to make sure it higher than elevation range. \n");
  }
  arr_out   = mxRealloc(arr_out,2*p*sizeof(int));   
  *num_out  = p;
}


/* ------ Get Difference Cells--------------------------------------------- 
  Get different cells of array 1 and array 2
  Getdiff will return the cells surround the old depression cells           
--------------------------------------------------------------------------*/
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
      arr_temp[K*p] 	= arr1[K*i];
      arr_temp[K*p+1]	= arr1[K*i+1];
      p++;
    }
  }
  arr_temp 	= realloc(arr_temp,2*p*sizeof(int));

  for (i=0;i<p-1;i++){
    for (j=i+1;j<p;j++){
      if (arr_temp[K*i] != -999 && arr_temp[K*i+1] != -999 && arr_temp[K*i]==arr_temp[K*j] && arr_temp[K*i+1]==arr_temp[K*j+1]){
        arr_temp[K*j] 	= -999;
        arr_temp[K*j+1] = -999;
      }
    }
  }
  for (i=0;i<p;i++){
    if (arr_temp[K*i] != -999 && arr_temp[K*i+1] != -999){
      arr_out[K*q] = arr_temp[K*i];
      arr_out[K*q+1] = arr_temp[K*i+1];
      q++;
    }
  }

  arr_out 	= mxRealloc(arr_out,2*q*sizeof(int));
  *num_out	= q;
}


/* ------ Get Neighbor Cells---------------------------------------------
	Get neighbor cells around point(s)
	D8 neighbor of all cellls
-------------------------------------------------------------------------*/
void getneigh(int arr_size, int arr[], int arr_out[])
{
  int p, q, K;

  K = 2;
  p = 0;
  arr_out = realloc(arr_out,16*arr_size*sizeof(int));

  for (q=0;q<arr_size;q++){
    arr_out[K*p+0]			= arr[K*q+0] - 1;
    arr_out[K*p+1]			= arr[K*q+1] - 1;
    arr_out[K*(p+1)+0]	= arr[K*q+0] - 1;
    arr_out[K*(p+1)+1]	= arr[K*q+1];

    arr_out[K*(p+2)+0]	= arr[K*q+0] - 1;
    arr_out[K*(p+2)+1]	= arr[K*q+1] + 1;
    arr_out[K*(p+3)+0]	= arr[K*q+0];
    arr_out[K*(p+3)+1]	= arr[K*q+1] - 1;

    arr_out[K*(p+4)+0]	= arr[K*q+0];
    arr_out[K*(p+4)+1]	= arr[K*q+1] + 1;
    arr_out[K*(p+5)+0]	= arr[K*q+0] + 1;
    arr_out[K*(p+5)+1]	= arr[K*q+1] - 1;

    arr_out[K*(p+6)+0]	= arr[K*q+0] + 1;
    arr_out[K*(p+6)+1]	= arr[K*q+1];
    arr_out[K*(p+7)+0]	= arr[K*q+0] + 1;
    arr_out[K*(p+7)+1]	= arr[K*q+1] + 1;
    p += 8;
  }
}



