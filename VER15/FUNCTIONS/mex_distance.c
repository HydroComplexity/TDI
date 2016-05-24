/*
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
  double *Cen_sing, *dist;
  double d;
  int i, j, k, low, up, num;

  /* Get value and dimensions of input matrix M */
  Cen_sing			= (double *)mxGetPr(prhs[0]); 														
  num					= mxGetM(prhs[0]); 
  k						= mxGetScalar(prhs[1]); 
  dist 				= 	(double*)malloc(num*sizeof(double));
  for (i=0;i<num;i++){
    dist[i] = 10000;
    if (i < k){
      low = i;
    }
    else{
      low = k;	
    }

    if (i > num-k){
      up = num-i;
    }
    else{
      up = k;
    }
    for (j=i-low;j<i+up;j++){
      if (i != j){
        d = pow(Cen_sing[0*num+i]-Cen_sing[0*num+j],2) + pow(Cen_sing[1*num+i]-Cen_sing[1*num+j],2);
        if (d < dist[i]){
          dist[i] = d;	
        }	
      }
    }
  }
  plhs[0] = mxCreateNumericMatrix(num,1, mxDOUBLE_CLASS, mxREAL);
  memmove(mxGetData(plhs[0]),dist,num*sizeof(double));

} /* end mexFunction() */
