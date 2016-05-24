function [A] = Pud_connect(D,n)

D(D>0)      = 1;
[r,t]       = size(D);
A           = zeros(r,t);
% Get connectivity of no-flow cells from image D.
[B,NUM]     = bwlabel(D,8);                                                 % value 8 represents double precision

[num_hist,med] = mex_regHistogram(B,NUM+1);                                 % fast histogram for connectivity image B

len_grenum  = length(find(num_hist >= n));                                  % identify numbers of connected objects have > 1 cells
len_lesnum  = length(find(num_hist < n));                                   % identify numbers of sing object has ONLY 1 cells

ind_hist    = find(num_hist >= n);

S           = sparse(B);                                                    % sparse matrix of B
linC        = unique(B);                                                    % Get values of matrix B without repetition.

mex_Pud_connect(S,A,linC,len_grenum,ind_hist); 