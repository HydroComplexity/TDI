function [Pud_pts,Puddle,ind_pud,ind_pudcum,Thres_pts,Threshold,ind_thres,ind_threscum,volume,area,max_depth,mean_depth]... 
          = pudsear(fid,MAT,MAT_ori,MAT_lev,Cen_multi,ind_mul,Cen_single,mink,max_size)
% 
% Find cells that belong to the depression and that are the thresholds.
% The searching process is started based on centers (or flats) found from 
%   the D8_center function.
% This function is inspired by the algorithm proposed by Chu et al (2011)
%-------------------------------------------------------------------------
%     Input arguments:
%       MAT:      a struct variable containing at least the fields:     
%         MAT.grid - a K x J matrix of elevations
%         MAT.dx   - grid spacing in x-direction
%         MAT.dy   - grid spacing in y-direction
%       MAT_ori:  Original DEM image in struct variable format.
%
%       Cen_pts   : Both single & multi center points arrays obtained 
%                     from D8_center
%       Cen_single: Single center points arrays obtained from D8_center
%       Cen_multi : Multi center points arrays obtained from D8_center
%       mink      : Minkowski threshold
%       max_size  : Maximum center size above which searching process will
%                     be implemented
%
%     Return arguments:
%       Pud_pts   :  Depression points arrays obtained from searching 
%       Puddle    :  a matrix of cells that are depression cells
%
%       Thres_pts :  Threshold  points arrays obtained from searching 
%       Threshold :  a matrix of cells that are threshold cells
%
%       volume    :  Array of corresponding volume for each puddle
%       area      :  Array of corresponding surface area for each puddle
%--------------------------------------------------------------------------
% PUDSEAR Function
% Created by: Phong Le
% Last Modified: 11-Nov-2012
%--------------------------------------------------------------------------

if (nargin ~= 9) 
  disp('ERROR: Please check the input variables.'); 
  disp('Inputs include the fid, DEM, DEM_ori, DEM_lev, 3 Center-point sets and mink');
  return 
end

Psi           = MAT.grid;                                                   % Elevation matrix
dx            = MAT.dx;                                                     % grid sizes
dy            = MAT.dy;
[M,N]         = size(Psi);                                                  % Matrix size
Pud_pts       = -999*ones(1,2);                                             % Initialize the Puddle points (1st column)
Thres_pts     = -999*ones(1,2);                                             % Initialize the Threshold points (1st column)
ind_pud       = -999*ones(1,1);                                             % Initialize the Threshold points (1st column)
ind_thres     = -999*ones(1,1);                                             % Initialize the Threshold points (1st column)
Puddle        = zeros(M,N);                                                 % Initialize binary Puddle matrix (0 & 1)
Threshold     = zeros(M,N);                                                 % Initialize binary Threshold matrix (0 & 1)

num_single    = size(Cen_single,1);                                         % Number of single center cells
num_multi     = length(ind_mul) -1;                                         % Number of multi center cells

%se           = [1 1; 1 0; 1 -1; 0 -1; 0 1; -1 -1; -1 0; -1 1];             % Structuring element for Minkowski operation (Dilation)

%% . . .DEPRESSION AND THRESHOLD IDENTIFICATION . . . . . . . . . . . . . .
%... MULTI-CENTER LOOP SEARCHING       
if ~isempty(Cen_multi)
  [Pud_mul, Thres_mul, ind_mul_pud, ind_mul_thres, Puddle, Threshold]... 
            = mulsear(Psi, Cen_multi, num_multi, ind_mul, ...
                      Puddle, Threshold, mink, max_size);                   % Call MULti SEARch matlab function
  %... Identify real centers after DEM filling for level > 1
  ind_mul_cum = cumsum(ind_mul);
  [Cen_real_mul] = mex_getmincenters(MAT_ori.grid,double(Cen_multi)-1,...
                                     ind_mul_cum);
  Center_real = Cen_real_mul;

  Pud_pts     = cat(1,Pud_pts,Pud_mul+1);                                   % Concatenate arrays
  Thres_pts   = cat(1,Thres_pts,Thres_mul+1);                               % Concatenate arrays
  ind_pud     = cat(1,ind_pud,ind_mul_pud(2:end));                          % Concatenate arrays
  ind_thres   = cat(1,ind_thres,ind_mul_thres(2:end));                      % Concatenate arrays
  
end

%... SINGLE-CENTER LOOP SEARCHING  
%... Identify & separate single centers & single thresholds
if ~isempty(Cen_single)
  Cen_pud_sing  = reshape(Cen_single',1,num_single*2)';                     
  max_elev      = 99999;%max(max(MAT_ori.grid));                            % Put maximum elevation to extreme high value
  [MAT_new]     = mex_fill_sing(Cen_single, MAT_lev.grid, max_elev);        % Fill-up single center cells to max elevation of DEM
  [TCA_sub, Bound_sub, Flow_dir_sub] ... 
      = mexD8(MAT_new, dy/dx, 0, -9999);                                   	% Get the new flow direction based on fill-up DEM

  [Cen_s, Thres_s, Cen_other, Puddle, Threshold]...
      = mex_sing_centhres(Cen_pud_sing, Flow_dir_sub, MAT_new,  ...
                          Puddle,Threshold);                                % Identify single threshold of single-center cells - No puddle except center.

  Cen_s         = Cen_s(Cen_s ~= 0);                                        % Remove 0-value at the end (shorter length)
  Cen_s         = reshape(Cen_s,2,length(Cen_s)/2)';
  Thres_s       = Thres_s(Thres_s ~= 0);                                  	% Remove 0-value at the end (shorter length)
  Thres_s       = reshape(Thres_s,2,length(Thres_s)/2)';  
  i_s           = 1:1:length(Cen_s);
  Cen_real_s    = zeros(length(Cen_s),1);  
  Cen_real_s    = MAT_ori.grid((Cen_s(i_s,1)-1)*N+Cen_s(i_s,2));
  Center_real   = cat(1,Center_real,Cen_real_s);
  
  Pud_s         = Cen_s;                                                    % Single-center, Pud cells are center cells
  Pud_pts       = cat(1,Pud_pts, Pud_s);                                    % Concatenate arrays
  ind_pud_add   = ones(size(Pud_s,1),1);
  ind_pud       = cat(1,ind_pud,ind_pud_add);                               % Concatenate arrays

  Thres_pts     = cat(1,Thres_pts,Thres_s);                                 % Combine multi & single centers & Thresholds lists
  ind_thres_add = ones(size(Thres_s,1),1);
  ind_thres     = cat(1,ind_thres,ind_thres_add);                           % Concatenate arrays

  Cen_other     = Cen_other(Cen_other ~= 0);                          
  num_other     = length(Cen_other)/2;
  Cen_other     = reshape(Cen_other,2,num_other)';
  ind_other     = ones(num_other+1,1);
  ind_other(1)  = 0;
  i_other       = 1:1:length(Cen_other);
  Cen_real_other= zeros(length(Cen_other),1);  
  Cen_real_other= MAT_ori.grid((Cen_other(i_other,1)-1)*N+Cen_other(i_other,2));
  Center_real   = cat(1,Center_real,Cen_real_other);

  %... Identify & separate single centers, puddle cells & single thresholds
  [Pudother_pts, Thresother_pts, ind_pudother, ind_thresother, ...
   Puddle, Threshold]... 
      = mulsear(Psi,    Cen_other,  num_other,  ind_other, ...
                Puddle, Threshold,  mink,       max_size);                  % Call MULti SEARch matlab function

  if length(ind_pudother) > 1            
    Pud_pts     = cat(1,Pud_pts, Pudother_pts+1);        
    Thres_pts   = cat(1,Thres_pts,Thresother_pts+1);                        % Combine multi & single centers & Thresholds lists
  
    ind_pud     = cat(1,ind_pud,ind_pudother(2:end));
    ind_thres   = cat(1,ind_thres,ind_thresother(2:end));  
  end
end
  
%% . . .ESTIMATE SURFACE AREA AND VOLUME. . . . . . . . . . . . . . . . . .
if ~isempty(size(Pud_pts,1)>=2)
  Pud_pts      	= Pud_pts(2:end,:);
  Thres_pts   	= Thres_pts(2:end,:);

  ind_pud(1)   	= 0;
  ind_pudcum   	= cumsum(double(ind_pud));
  ind_thres(1)	= 0;
  ind_threscum  = cumsum(double(ind_thres));
  area          = ones(length(ind_pud)-1,1);
  area          = ind_pud(2:end)*dx*dy;
  
  [volume,max_depth] ... 
      = mex_geometric(MAT_ori.grid, double(Pud_pts)-1, double(Thres_pts)-1,...
                      ind_pudcum, ind_threscum, Center_real, dx, dy);
  mean_depth   	= volume./double(area);
else
  area         	= nan;
  volume        = nan;
  max_depth   	= nan;
  mean_depth  	= nan;
end

%% . . .STATISTICS AND PRINTING . . . . . . . . . . . . . . . . . . . . . .
disp(       ['       . . .DEPRESSION & THRESHOLD STATISTICS . . . . .'                  ]);
disp(       ['       Number of depressions: ',num2str(length(ind_pud)-1)                ]);
disp(       ['       Number of cells in the largest depression: ',num2str(max(ind_pud)) ]);
disp(       ['       Number of thresholds: ',num2str(length(ind_pud)-1)                 ]);
disp(       ['       Largest # of thresholds in one puddle: ',num2str(max(ind_thres))   ]);

fprintf(fid,['      . . .DEPRESSION & THRESHOLD STATISTICS . . . . .\n'                 ]);
fprintf(fid,['      Number of depressions: ',num2str(length(ind_pud)-1),'\n'            ]);
fprintf(fid,['      Number of cells in the largest depression: ',num2str(max(ind_pud)),'\n']);
fprintf(fid,['      Number of thresholds: ',num2str(length(ind_pud)-1),'\n'             ]);
fprintf(fid,['      Largest # of thresholds in one puddle: ',num2str(max(ind_thres)),'\n']);
