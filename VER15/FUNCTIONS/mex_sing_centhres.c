/*=======================================================================
* mex_sing_centhres.c - C MEX file used to identify the single-theshold
*     of single-center puddles.
*------------------------------------------------------------------------
*     usage: [Cens,Thres,Pud_out,Thres_out] 
              = mex_fill_sing(Center,M,maxE,Pud_in,Thres_in)
*------------------------------------------------------------------------
*     Input arguments:
*       Cen:        list (x & y) locations of single-center cells
*       M:          a K x J matrix of flow direction
*       DEM:        a K x J matrix of elevation
*
*     Return arguments:
*       Cens:       list (x & y) locations of single-center cells
*       Thres:      list (x & y) locations of single-threshold cells
*       Pud_out:    a K x J matrix of puddle cells
*       Thres_out:  a K x J matrix of threshold cells
*------------------------------------------------------------------------
* This is a MEX-file for MATLAB.
* Created by: Phong Le <phongle1@illinois.edu>
* This code may be freely redistributed in its entirety provided that 
  this copyright notice is	not removed.
*========================================================================
*/
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize n;
  double *M, *DEM, *A, *Cens, *Cens_other, *Thres, *Pud_in, *Pud_out, *Thres_in, *Thres_out;			
  int i, j, p, q, K, J, row, col, row_min, col_min, colup, coldown;
  int num_loop, zero_flag, ind[8][2]; 
  double *Cen, min_array;
  int Outdims[2];               /* dimensions of the output arrays */

  /* Get pointers to inputs, assign default values if not supplied */
  n   = mxGetNumberOfElements(prhs[0]);
  Cen = mxGetPr(prhs[0]);

  /* Get value and dimensions of input matrix M */
  M   = (double *)mxGetPr(prhs[1]); 
  K   = mxGetM(prhs[1]);
  J   = mxGetN(prhs[1]);	

  /* Get DEM value for other inputs */
  DEM = (double *)mxGetPr(prhs[2]); 				

  Outdims[0]  = 1; 
  Outdims[1]  = n; 
  num_loop    = n/2;

  /* Create arrays for the return arguments*/
  Cens        = mxGetPr(plhs[0] = mxCreateNumericArray(2, Outdims, mxDOUBLE_CLASS, mxREAL));
  Thres       = mxGetPr(plhs[1] = mxCreateNumericArray(2, Outdims, mxDOUBLE_CLASS, mxREAL));
  Cens_other  = mxGetPr(plhs[2] = mxCreateNumericArray(2, Outdims, mxDOUBLE_CLASS, mxREAL));
  Pud_out     = (double *)mxGetPr(plhs[3] = mxDuplicateArray(prhs[3]));
  Thres_out   = (double *)mxGetPr(plhs[4] = mxDuplicateArray(prhs[4]));
  
  /* Update max elevation to matrix M */
  p = 0;
  q = 0;
  for (i=0; i<num_loop; i++) {
    col       = (int) Cen[2*i]-1;
    row       = (int) Cen[2*i+1]-1;
    colup     = K * (col+1);
    coldown   = K * (col-1);

    ind[0][0] = col - 1;
    ind[0][1] = row - 1;
    ind[1][0] = col - 1;
    ind[1][1] = row;
    ind[2][0] = col - 1;
    ind[2][1] = row + 1;
    ind[3][0] = col;
    ind[3][1] = row - 1;
    ind[4][0] = col;
    ind[4][1] = row + 1;
    ind[5][0] = col + 1;
    ind[5][1] = row - 1;
    ind[6][0] = col + 1;
    ind[6][1] = row;
    ind[7][0] = col + 1;
    ind[7][1] = row + 1;

    zero_flag = 1;
    for (j=0;j<8;j++){
      if (M[K*ind[j][0] + ind[j][1]] == 0){
        zero_flag = 0;
        break;
      }
    }

    if (zero_flag == 1){
      min_array = DEM[K*ind[0][0]+ind[0][1]];
      col_min   = ind[0][0];
      row_min   = ind[0][1];
      for (j=1; j<8;j++) {
        if (DEM[K*ind[j][0]+ind[j][1]] < min_array){
          min_array = DEM[K*ind[j][0]+ind[j][1]];
          col_min   = ind[j][0];
          row_min   = ind[j][1];
        }
      }
      Thres[2*p]    = col_min + 1;
      Thres[2*p+1]  = row_min + 1;
      Thres_out[K*col_min+row_min] = 1;
      
      Cens[2*p]     = col + 1;
      Cens[2*p+1]   = row + 1;
      Pud_out[K*col+row] = 1;	
      p++;
    }
    else {
      Cens_other[2*q]   = col + 1;
      Cens_other[2*q+1] = row + 1;
      q++;
    }
  }
} /* end mexFunction() */
