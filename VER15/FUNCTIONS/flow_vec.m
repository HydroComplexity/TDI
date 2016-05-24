function [DX,DY] = flow_vec(MAT,flowDir)

Psi   = MAT.grid;
[M,N] = size(Psi);
DX    = zeros(M,N);
DY    = zeros(M,N);

for i = 1:M
  for j =1:N
%    indx = flowDir{i,j}(2);
%    indy = flowDir{i,j}(1);    
    DX(i,j) = (flowDir{i,j}(2) - j)-0.5;
    DY(i,j) = (flowDir{i,j}(1) - i)-0.5;
  end
end