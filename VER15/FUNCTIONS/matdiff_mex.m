function [C,D] = matdiff(A,B)
% MATRIX DIFFERENCE: Remove matrix B from matrix A if similar values are
% found between two matrices
%
%   Given matrix A with size = M x 2, and matrix B with size = N x 2
%   The matdiff will remove any values in A that are found in B (both
%   columns 1 & 2 in A must be matched with columns in B). So far, the
%   function only works for matrices with 2 columns.
%
%-------------------------------------------------------------------------
%   Example:
%   Inputs:
%     A = [1 2; 3 1; 2 3; 4 5; 3 6; 7 2]; 
%
%     B = [3 1; 4 5; 6 9];
%
%   Output:
%     C = [[1 2; 2 3; 3 6; 7 2]; 
%
%-------------------------------------------------------------------------
% MATDIFF
% Created by: Phong Le
% University of Illinois at Urbana-Champaign
% Email: phongle1@illionis.edu
% Last Modified: 09-Nov-2012
%--------------------------------------------------------------------------
 
if nargin<2
    error('Error: Requires two inputs');
end

[M,N ] = size(B);
B = reshape(B,N,M)';

szA=size(A); szB=size(B);
if szA(2)~= 2 || szB(2)~=2 || numel(szA)>2 || numel(szB)>2
    error('Input must be Mx2 and Nx2 arrays');
end

ind = ismember(A,B,'rows');
A = A(ind==0,:);

Anew            = unique(A,'rows');  
nanRows         = any(isnan(Anew), 2);
Anew(nanRows,:) = [];
C = int32(Anew');
D = int32(size(C,2));