function [B, NUM, A_multi, ind_mul, A_single, ind_single] = D8_checkzero(D,DEM0,mink)
%
% Identifies the lowest points of depressions that may become puddles and
% computes upslope contributing area.
%-------------------------------------------------------------------------
%     Input arguments:
%       D:      	a matrix includes flow direction
%       DEM0:     a matrix includes 0 values from mexD8.
%       mink:     Minkowski threshold
%
%     Return arguments:
%       A:        a list includes sets of center(s)' location. The first
%                   column is unused. Every next two columns represent the
%                   position (x,y) of cells belong to a puddle.
%                   Denote m is the width of matrix A
%                   The total number of depressions is: (m-1)/2;                 
%
%       B:        a matrix of connectivity of values 0
%
%       NUM:      Number of connected objects in image DEM0
%
%       A_multi:  a list includes sets of multi (connected) centers'
%                   location. It is a part of output A
%
%       A_single: a list includes sets of single center locations. It is
%                   also a part of output A. (A = A_multi + A_single)
%--------------------------------------------------------------------------
% FUNCTION D8_checkzeo
% Created by: Phong Le
% Last Modified: 10-Nov-2012
%--------------------------------------------------------------------------

%... Check the number of inputs
if (nargin ~= 3) 
  disp('Check the inputs. . . . .');
  disp('The flow direction matrix , DEM, and mink_thres must be included');
  return 
end

[M,N]         = size(D);                                                    % Get size of matrix
B             = zeros(M,N);                                                 % Initializing zero matrix for faster bwlabel

%... Flow direction is from 1-8, CCW starting from the East. Value 0 means no flow. 
%... The commands below convert non-zero values to 0, but 0 to non-zero values (9)
D(D==0)       = 9;
D(D<9)        = 0;

%... Get connectivity of no-flow cells from image D.
[B,NUM]       = bwlabel(D,8);                                               % value 8 represents double precision
Bs            = B;
DEM_one       = reshape(DEM0, M*N,1);                                       % convert DEM0 into single column
A_multi       = -999*ones(1,2);                                             % initializing output A (two-column matrix)

%... Structuring element for dilation process
se            = [1 1; 1 0; 1 -1; 0 -1; 0 1; -1 -1; -1 0; -1 1];

[num_hist,med]= mex_regHistogram(B,NUM+1);                                  % Fast C histogram for connectivity image B

len_grenum    = length(find(num_hist > 1));                                	% identify numbers of connected objects that have > 1 cells
len_lesnum    = length(find(num_hist == 1));                                % identify numbers of sing object has ONLY 1 cells
ind_hist      = find(num_hist > 1);

S             = sparse(B);                                                  % sparse matrix of B
linC          = unique(B);                                                  % Get values of matrix B without repetition.
ind_mul       = zeros(len_grenum,1);                                        % Initialize index of multi center
ind_mul(1)    = 1;    

%%. . .LOOP FOR SEARCHING MULTI-CENTERS.
[A_multi]     = mex_mul_centhres(DEM0,S,linC,len_grenum,ind_hist,...
                                 ind_mul,Bs,mink);
A_multi       = double(A_multi)'+1;                                         % Add 1 to all value for converting back from C to MATLAB
A_multi       = A_multi(2:end,:);
ind_mul       = ind_mul(ind_mul ~= 0);                                      % Remove index_multi (represent # of cells in flat) equal to 0
ind_mul(1)    = 0;                                                          % First index equal to 0 for cumulative later

%%. . .LOOP FOR SEARCHING SINGLE-CENTERS.
[A_add]       = mex_single(Bs,len_lesnum);                                  % Get locations of single center - edges are not included.
A_non0        = A_add(A_add~=0);                                            % Remove zero values of A_add (single).
A_single      = reshape(A_non0,2,length(A_non0)/2)';                        % Reshape to two-column matrix
ind_single    = ones(length(A_single)+1,1);                                 % index = 1 (because single) for all center.
ind_single(1) = 0;                                                          % First index equal to 0 for cumulative later

