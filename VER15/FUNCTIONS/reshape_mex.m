function [B] = reshape_mex(A)
  [M,N ] = size(A);
  B = reshape(A,N,M)';
end
