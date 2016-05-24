function [MAT_cor] = DEM_cor(MAT,Cen_multi,ind_mul,Cen_single,nodata,mink)
%
% [Pud_pts, Puddle, Thres_pts, Threshold, volume, area, Cens, Thres] 
%         = pudsear(MAT, MAT_ori, Cen_pts, Cen_single, Cen_multi)
% 
%-------------------------------------------------------------------------
%     Input arguments:
%       MAT:      a struct variable containing at least the fields:     
%         MAT.grid - a K x J matrix of elevations
%         MAT.dx   - grid spacing in x-direction
%         MAT.dy   - grid spacing in y-direction
%       MAT_ori:  Original DEM image in struct variable format.
%
%--------------------------------------------------------------------------
% PUDSEAR Function
% Created by: Phong Le
% Last Modified: 11-Nov-2012
%--------------------------------------------------------------------------

if (nargin ~= 6) 
  disp('ERROR: Please check the input variables.'); 
  disp('Inputs include the DEM, DEM_ori, DEM_lev and 3 Center-point sets');
  return 
end

Psi           = MAT.grid;                                                     % Elevation matrix
dx            = MAT.dx;                                                       % grid sizes
dy            = MAT.dy;

%... MULTI-CENTER CORRECTION
if ~isempty(Cen_multi)
  ind_cum = cumsum(ind_mul);
  [DEM_correct_mul] = mex_DEM_correct_mul(Psi,Cen_multi'-1,ind_mul,ind_cum,nodata,mink);
  MAT_cor.grid      = DEM_correct_mul;  
  Psi               = DEM_correct_mul;  
end
  
%... SINGLE-CENTER CORRECTION
if ~isempty(Cen_single)  
  [DEM_correct_sing] = mex_DEM_correct_sing(Psi,Cen_single'-1,nodata);
  MAT_cor.grid       = DEM_correct_sing; 
end

MAT_cor.dx = dx;
MAT_cor.dy = dy;
