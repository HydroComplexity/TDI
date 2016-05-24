#include "mex.h"
#include <string.h>
#include "matrix.h"
#include <math.h>

void getminthres(double MAT[], int M, int N, int arr_in[], int size_in, double *min_out)
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

void getneighthres(int arr[], int arr_out[])
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *DEM, *Threshold, *Pud_pts, *ind_pudcum, minval;
  int *At, *Dt;
  int i, j, M, N, Mp, num_pud;

  /* Get value and dimensions of input matrix M */
  DEM       = (double *)mxGetPr(prhs[0]);
  M         = mxGetM(prhs[0]);
  N         = mxGetN(prhs[0]);	
  Threshold = (double *)mxGetPr(prhs[1]);

  Dt = (int*)malloc(16*sizeof(int));
  At = (int*)malloc(2*sizeof(int));

  for (i=1;i<M-1;i++){
    for (j=1;j<N-1;j++){
      if (Threshold[M*j+i] > 1){
        At[0] = j;
        At[1] = i;
        getneighthres(At,Dt);
        getminthres(DEM, M, N, Dt, 8, &minval);
        if (DEM[M*j+i] > minval)
          Threshold[M*j+i] = 1;
      }
    }
  }

} /* end mexFunction() */
