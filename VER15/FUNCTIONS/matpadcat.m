function [C] = matpadcat(A,B)
% MATPADCAT - concatenate multi-column matrices with different sizes by 
%             padding with NaN.
%
%   M = MATPADCAT(M1, M2) concatenates the matrix M1 and M2 into one large
%   matrix. The matrices do not need to have the same size. More columns 
%   will be added to M1 and combine with M2 as illustrated below.
%   For padding by rows, use matrix transpose for your purposes.
%
%-------------------------------------------------------------------------
%   Example:
%   Inputs:
%     A = [1 1 1; 2 2 2; 3 3 3];
%     B = [4 4; 5 5];
%
%   Output:
%     C = 1   1   1   4   4
%         2   2   2   5   5
%         3   3   3   NaN NaN
%   
%--------------------------------------------------------------------------
%   This work is inspired by the PADCAD function (See information below for 
%   more details on PADCAT).
%     PADCAT: version 1.2 (oct 2011)
%     (c) Jos van der Geest
%     email: jos@jasen.nl
%--------------------------------------------------------------------------
% MATPADCAT
% Created by: Phong Le
% University of Illinois at Urbana-Champaign
% email: phongle1@illionis.edu
% Last Modified: 10-Nov-2012
%--------------------------------------------------------------------------

if nargin<2
    error('Requires two inputs');
end

szA   = size(A);    szB   = size(B);
legA  = szA(1);     widA  = szA(2); 
legB  = szB(1);     widB  = szB(2); 

dif_size = abs(legA - legB);

if legA > legB
  Badd  = nan(dif_size,widB);
  Bnew  = cat(1,B,Badd);
  C     = cat(2,A,Bnew);
elseif legA < legB
  Aadd = nan(dif_size,widA);
  Anew  = cat(1,A,Aadd);
  C     = cat(2,Anew,B);
else
  C     = cat(2,A,B);
end