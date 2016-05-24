function [MAT_new] = threshold_fill(MAT,MAT_pud,Pud_pts,ind_pudcum,Thres_pts,ind_threscum,Threshold)
%
% [DEM_new] = threshold_fill(DEM,Cen_pts, Thres_pts)
% 
% Search shared common threshold among depression. If shared thresholds are
% found, fill the related depression to the threshold elevation.
% Depression and Threshold points from D8_Center and pudsear functions are
% used for filling.
% 
%-------------------------------------------------------------------------
%     Input arguments:
%       MAT:      a struct variable containing at least the fields
%         MAT.grid: a K x J matrix of elevations
%         MAT.dx:   grid spacing in x-direction
%         MAT.dy:   grid spacing in y-direction
%
%       Pud_pts   :  Depression points arrays obtained from searching 
%       Thres_pts :  Threshold  points arrays obtained from searching 
%       Thresdhol :  Matrix of Threshold points
%
%     Return arguments:
%       MAT_new:    a new struct variable containing at least the fields
%         MAT_new.grid: a K x J matrix of filling threshold elevations
%         MAT_new.dx:   grid spacing in x-direction
%         MAT_new.dy:   grid spacing in y-direction
%
%--------------------------------------------------------------------------
% DEM_FILL
% Created by: Phong Le
% University of Illinois at Urbana-Champaign
% email: phongle1@illionis.edu
% Last Modified: 14-Nov-2012
%--------------------------------------------------------------------------

if (nargin ~= 7) 
  disp('Check the DEM inputs. Puddle and Threshold points are required');
  return 
end
Psi         = MAT.grid;                                                     % Elevation matrix

%... Check if threshold is eligible for merging
mex_check_thres(MAT_pud,Threshold);

[i,j]       = find(Threshold > 1);
C           = unique([j,i],'rows');
D           = ismember(Thres_pts,C,'rows');
ind_one     = find(D==1);


if ~isempty(ind_one)
  mex_fill_threshold(Psi, double(Pud_pts)-1, double(Thres_pts)-1, ind_pudcum, ind_threscum, double(ind_one)-1);
else
  disp(['Empty ind_one']);
end

MAT_new.grid  = Psi;
MAT_new.dx    = MAT.dx;
MAT_new.dy    = MAT.dy;












