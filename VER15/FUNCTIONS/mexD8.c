/*
* ========================================================================
* mexD8.c - C MEX file to compute upslope contributing area using the D8
*           (steepest descent) algorithm of O'Callaghan & Mark (1984)
*
*     usage: [A, B, D] = mexD8(M,dy/dx,flood,ndv)
*
*     Input arguments:
*       M      a K x J matrix of elevations
*       dy/dx  the ratio of grid spacings in the y and x directions. e.g.,
*              if dy=4m and dx=2m, dy/dx=2. (optional: default = 1)
*       flood  if 1, routes flow through local minima in M (potentially
*              time-consuming). if zero, skips this step (faster).
*              (optional: default = 0)
*       ndv    matrix value that indicates no data. this will be used to
*              identify the boundaries of your elevation matrix. If all
*              the elements in your matrix have valid elevations, just set
*              ndv to some value that does not appear in your elevations
*              (e.g., -1 or 9999). (optional: default = 0)
*
*     Return arguments:
*       A      a matrix of total contributing areas (in units of CELLS)
*       B      a matrix of cells that receive drainage from boundaries
*       D      a matrix of D8 drainage directions (1-8, CCW from east, zero
*              at sinks and boundaries) 
*
* This is a MEX-file for MATLAB.
*
* Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
* 
* This program is free software: you can redistribute it and/or modify it 
* under the terms of the GNU General Public License as published by the 
* Free Software Foundation. You should have received a copy of the GNU 
* General Public License along with this program.  If not, see 
* http://www.gnu.org/licenses.
* ========================================================================
*/

#include "mex.h"
#include "matrix.h"
#include <math.h>
  
void GetBound(const int K, const int J, double M[], double Bdy[], double B[], double *ndv)
{

int i, j, Kj, Kju, Kjd;
double n;

n = *ndv;

for (j=0; j<J; j++) {
    Kj=K*j;
    Kju=K*(j+1);
    Kjd=K*(j-1);
    for (i=0; i<K; i++) {
        if (M[Kj+i] != n) { /* if M(i,j) has data */
            if (i==0 || i==K-1 || j==0 || j==J-1) { /* if on matrix boundary */
                Bdy[Kj+i] = 2;
                B[Kj+i] = 1; /* mark as boundary-influenced */
            } else if (M[Kju+i] == n || M[Kju+i-1] == n || M[Kj+i-1] == n || M[Kjd+i-1] == n || M[Kjd+i] == n || M[Kjd+i+1] == n || M[Kj+i+1] == n || M[Kju+i+1] == n) { /* if any of its neighbors has nodata */
                Bdy[Kj+i] = 2;
                B[Kj+i] = 1; /* mark as boundary-influenced */
            } else { /* otherwise it's an interior elevation cell */
                Bdy[Kj+i] = 1;
            }
        } /* if M(i,j)==n, leave Bdy=0 */
    }
}
    
} /* end GetBound() */


void D8Dir(const int K, const int J, double M[], double Bdy[], double D[], double minima[], double *a, int *numMin)
{
int i, j, k, Kj, Kjup, Kjdown, theD;
double invdx, invdy, invdiag, theS, maxS;
double e0;

/* For boundaries and corners, restrict facets to the following: */
int numfacets, *facetidx;
int ULfacets[]	={6,7,0};
int URfacets[]	={4,5,6};
int LLfacets[]	={0,1,2};
int LRfacets[]	={2,3,4};
int Lfacets[]	={6,7,0,1,2};
int Rfacets[]	={2,3,4,5,6};
int Tfacets[]	={4,5,6,7,0};
int Bfacets[]	={0,1,2,3,4};
int Allfacets[]	={0,1,2,3,4,5,6,7};

invdx	= 1;
invdy	= 1;
invdiag	= 1;

  for (j=0; j<J; j++) { /* Loop through columns */

      Kj=K*j;
      Kjup=K*(j+1);
      Kjdown=K*(j-1);

      for (i=0; i<K; i++) { /* Loop through rows */
          
          if (Bdy[Kj+i] > 0) { /* if it is nodata, skip it */

              e0=M[Kj+i]; /* Elevation at (i,j) */
              maxS = 0; /* Set max slope to zero */
              theD = 0;

              if (i>0 && i<(K-1) && j>0 && j<(J-1)) { /* Interior cell */
                  numfacets=8;
                  facetidx=Allfacets;
              } else if (i==0 && j==0) { /* UL corner */
                  numfacets=3;
                  facetidx=ULfacets;
              } else if (i==0 && j==J-1) { /* UR corner */
                  numfacets=3;
                  facetidx=URfacets;
              } else if (i==K-1 && j==0) { /* LL corner */
                  numfacets=3;
                  facetidx=LLfacets;
              } else if (i==K-1 && j==J-1) { /* LR corner */
                  numfacets=3;
                  facetidx=LRfacets;
              } else if (j==0) { /* Left boundary */
                  numfacets=5;
                  facetidx=Lfacets;
              } else if (j==J-1) { /* Right boundary */
                  numfacets=5;
                  facetidx=Rfacets;
              } else if (i==0) { /* Top boundary */
                  numfacets=5;
                  facetidx=Tfacets;
              } else if (i==K-1) { /* Bottom boundary */
                  numfacets=5;
                  facetidx=Bfacets;
              }

              for (k=0; k<numfacets; k++) { /* Loop through facets, ignoring any that have nodata */
                  switch (facetidx[k]) {
                      case 0:
                        theS=e0-M[Kjup+i];
                        if (Bdy[Kjup+i] == 0) {
                            theS = -1;
                        }
                        break;
                      case 1:
                        theS=(e0-M[Kjup+(i-1)])*invdiag;
                        if (Bdy[Kjup+(i-1)] == 0) {
                            theS = -1;
                        }
                        break;
                      case 2:
                        theS=(e0-M[Kj+(i-1)])*invdy;
                        if (Bdy[Kj+(i-1)] == 0) {
                            theS = -1;
                        }
                        break;
                      case 3:
                        theS=(e0-M[Kjdown+(i-1)])*invdiag;
                        if (Bdy[Kjdown+(i-1)] == 0) {
                            theS = -1;
                        }
                        break;
                      case 4:
                        theS=e0-M[Kjdown+i];
                        if (Bdy[Kjdown+i] == 0) {
                            theS = -1;
                        }
                        break;
                      case 5:
                        theS=(e0-M[Kjdown+(i+1)])*invdiag;
                        if (Bdy[Kjdown+(i+1)] == 0) {
                            theS = -1;
                        }
                        break;
                      case 6:
                        theS=(e0-M[Kj+(i+1)])*invdy;
                        if (Bdy[Kj+(i+1)] == 0) {
                            theS = -1;
                        }
                        break;
                      case 7:
                        theS=(e0-M[Kjup+(i+1)])*invdiag;
                        if (Bdy[Kjup+(i+1)] == 0) {
                            theS = -1;
                        }
                        break;
                  }


                  /* If this is a downslope facet AND the steepest yet, record its slope and direction */
                  if (theS>maxS) {
                      maxS=theS;
                      theD=facetidx[k]+1;
                  }

              } /* end looping through facets */

              /* If (i,j) is a local minimum (maxS remains zero), it retains a D8 direction of zero and is flagged in minima */
              if (maxS==0) {
                  minima[Kj+i]=1;
                  (*numMin)++;
              }
              D[Kj+i]=theD; /* assign D8 direction */
          }
      } 
  } 
} /* end D8Dir() */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  mxArray *Rptr, *Optr;
  double *M, *Bdy, *a, *fl, *ndv, *A, *B, *D, *F, *R, *O, *minima;
  int i, j, K, J, *numU, *numMin, outletidx[2], ndims=3, dims[]={0,0,8};

  /* Argument checking */
  if (nrhs < 1 || nrhs > 4) {
    mexErrMsgTxt("Wrong number of input arguments. Usage: [area boundary] = mexD8(elev,dy/dx,doflooding,nodataval)");
  } else if (nlhs > 3) {
    mexErrMsgTxt("Too many output arguments. Usage: [area boundary] = mexD8(elev,dy/dx,doflooding,nodataval)");
  }

  if (mxGetNumberOfDimensions(prhs[0]) != 2) {
    mexErrMsgTxt("Elevation array must be 2D.");
  }
  
  /* Get pointers to inputs, assign default values if not supplied */
  M = (double *)mxGetPr(prhs[0]); /* elevations */


  if (nrhs < 2) { 
      a = (double *)mxGetPr(mxCreateDoubleMatrix(1, 1, mxREAL));
      *a = 1; /* default */
  } else {
      a = (double *)mxGetPr(prhs[1]); /* dy/dx */
  }
  
  if (nrhs < 3) { 
      fl = (double *)mxGetPr(mxCreateDoubleMatrix(1, 1, mxREAL));
      *fl = 0; /* default */
  } else {
      fl = (double *)mxGetPr(prhs[2]); /* flood local minima, or not */
  }

  if (nrhs < 4) { 
      ndv = (double *)mxGetPr(mxCreateDoubleMatrix(1, 1, mxREAL));
      *ndv = 0; /* default */
  } else {
      ndv = (double *)mxGetPr(prhs[3]); /* nodata value */
  }

  /* Get dimensions of input matrix of elevations */
  K=dims[0]=mxGetM(prhs[0]);
  J=dims[1]=mxGetN(prhs[0]);
  
  /* Create arrays for the return arguments */
  A = (double *)mxGetPr(plhs[0]= mxCreateDoubleMatrix(K, J, mxREAL)); /* matrix of Total Contributing Areas*/
  B = (double *)mxGetPr(plhs[1]= mxCreateDoubleMatrix(K, J, mxREAL)); /* indicates which cells receive drainage from one or more boundary cells*/
  D = (double *)mxGetPr(plhs[2]= mxCreateDoubleMatrix(K, J, mxREAL)); /* matrix of D8 drainage directions*/

  /* Create internally used arrays */
  Bdy 		= (double *)mxGetPr(mxCreateDoubleMatrix(K, J, mxREAL)); /* matrix of boundary flags */
  minima 	= (double *)mxGetPr(mxCreateDoubleMatrix(K, J, mxREAL)); /* 1 if an element is a local minimum, zero otherwise */
  F 			= (double *)mxGetPr(mxCreateDoubleMatrix(K, J, mxREAL)); /* matrix that indicates what was flooded */
  numMin 	= mxCalloc(1, sizeof(int)); /* number of local minima */

  /* Find boundary elements, defined as cells with valid elevations that
   	are on the edge of the grid or have at least one neighbor with nodata. 
   	Bdy == 1 for interior points, 2 on boundaries, zero for nodata */
  GetBound(K,J,M,Bdy,B,ndv);

  /* Calculate D8 drainage directions */
  *numMin = 0;

  D8Dir(K,J,M,Bdy,D,minima,a,numMin); /* Calculate D8 drainage directions */
} /* end mexFunction() */
