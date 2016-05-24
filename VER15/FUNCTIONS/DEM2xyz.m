function DEM2xyz(MAT,filename)

if nargin~=2
    error('Error: Requires two inputs');
end

Psi         = MAT.grid;
[M,N]       = size(Psi);                                                    % Matrix size

fileID = fopen(filename,'w');

for i = 1:N    
    for j = 1:M
        x = i;
        y = j;
        z = Psi(j,i);
        fprintf(fileID,'%d %d %3.1f\n',x,y,z);        
    end
end

fclose(fileID);