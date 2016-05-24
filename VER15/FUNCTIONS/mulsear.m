function [Pud_pts, Thres_pts, ind_pud, ind_thres, Puddle, Threshold]... 
          = mulsear(MAT, Cen_pts, num_pts, ind_pts, Pud_in, Thres_in, mink,max_size)
%
% [Pud_pts, Thres_pts, Puddle, Threshold]
%         = mulsear(MAT_ori, Cen_pts, Pud_in, Thres_in)
% 
% Search puddle cells of non-single puddle&center type
%-------------------------------------------------------------------------
%     Input arguments:
%       MAT_ori:  Original DEM image with elevation, containing at 
%                 least the fields:     
%         .grid   - a K x J matrix of elevations
%         .dx     - grid spacing in x-direction
%         .dy     - grid spacing in y-direction
%
%       Cen_pts   : Both single & multi center points arrays obtained 
%                     from D8_center
%       Pud_in    : Image of Puddle cells from previous search
%       Thres_in  : image of Threshold cells from previous search
%       mink      : Minkowski threshold
%
%
%     Return arguments:
%       Pud_pts   :  Depression points arrays obtained from searching 
%       Puddle    :  a matrix of cells that are depression cells
%
%       Thres_pts :  Threshold  points arrays obtained from searching 
%       Threshold :  a matrix of cells that are threshold cells
%--------------------------------------------------------------------------
% MULSEAR Function
% Created by: Phong Le
% Last Modified: 02-Feb-2012
%--------------------------------------------------------------------------
if (nargin ~= 8) 
  disp('ERROR: Please check the input variables.'); 
  return 
end

Puddle    = Pud_in;
Threshold = Thres_in;
ind_cum   = cumsum(ind_pts);

[Pud_pts, Thres_pts,ind_pud, ind_thres]=mex_others_run(MAT,Cen_pts'-1,num_pts,ind_pts,ind_cum,Puddle,Threshold,mink,max_size);

Pud_pts   = reshape(Pud_pts,2,length(Pud_pts(:))/2)';
Thres_pts = reshape(Thres_pts,2,length(Thres_pts(:))/2)';