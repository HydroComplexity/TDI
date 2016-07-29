%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%%                           MAIN PROGRAM                                %%
%%             Depression & Puddle Identification Model                  %%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
% This model is used to identify                                          % 
%   + Depression locations & sizes (Center, Depression & Threshold cells) %
%   + Flow Direction                                                      % 
%   + Accumulative Area                                                   %
%-------------------------------------------------------------------------%
%   Created by  : Phong Le                                                %
%   Date        : November 12, 2012                                       %
%-------------------------------------------------------------------------%
%                                                                         
  clear all                         % Clear all memory and variables      
  close all                         % Close all programs and functions    
  clc
%
%% . . . . . . . . . . . . MODEL STRUCTURE. . . . . . . . . . . . . . . . .
%
% SWITCHES    :   Model conditional switches                              
% CONSTANTS   :   Indepentdent constants, unit conversions                
% PARAMETERS  :   Model Parameters                                        
%
%-------------------------------------------------------------------------%
%*************************************************************************%
%%                          USER SPECIFICATIONS                          %%
%   
  SWITCHES.save_on      = 1;      % 1 = Saving variables to .mat file, 
%                                 % 0 = No saving
%
  SWITCHES.plots_on     = 0;      % 1 = Plotting results, 
%                                 % 0 = No plotting
%
  SWITCHES.surfer       = 0;      % 1 = Export to Surfer GRD file, 
%                                 % 0 = No exporting
%
  SWITCHES.block        = 0;      % 1 = Split DEM into blocks quick search, 
%                                 % 0 = No block applied
%
  SWITCHES.correct      = 0;      % 1 = Apply correction for DEM, 
%                                 % 0 = No correction applied
%
  interpolate_method    = 0;      % 1 = Get nearest neighbor value 
                                  % 0 = inpainting interpolation, degree 3
%                                  
%%                      END OF USER SPECIFICATIONS                       %%
%*************************************************************************%
%
%% . . .ADD LIBRARY FUNCTION PATHS. . . . . . . . . . . . . . . . . . . . .
  addpath('./FUNCTIONS/');
  %...Check if RESULTS folder exits
  if (exist('RESULTS','dir') ~= 7)
    mkdir('RESULTS');
  end
  
%% . . . SETTING UP PARAMETERS. . . . . . . . . . . . . . . . . . . . . .
  ndv           = -9999;
  flood         = 0;  
  mink          = 10000;
  max_size      = 20000;  
  
% 
%% . . .LOAD TOPOGRAPHIC DATA . . . . . . . . . . . . . . . . . . . . . . .  
  file_in       = 'BPNM_DEM2';
  ex_file       = 'tif';
  %full_file_in  = ['../LiDAR_Data/BPNM/',file_in,'.',ex_file];
  full_file_in  = ['../LiDAR_Data/BPNM/BPNM_DEM2.tif'];
  [Z, R, bbox]  = geotiffread(full_file_in);
  info = geotiffinfo(full_file_in);
  [My,Mx]       = size(Z);        
  Z_in          = double(Z);
  [My_in,Mx_in] = size(Z_in);            
  
  if (SWITCHES.block)
    block_row = 1;
    block_col = 1;
    num_run   = [1,1];    
    if (num_run(1)>block_row || num_run(2)>block_col)
      disp(['Error: Block positions should be inside the domain']);
      return
    end
    len_row   = floor(My/block_row);
    len_col   = floor(Mx/block_col);
    DEM_in    = double(Z_in((num_run(1)-1)*len_row+1:num_run(1)*len_row,...
                            (num_run(2)-1)*len_col+1:num_run(2)*len_col));
    file_out  = [file_in,'_',num2str(num_run(1)),'_',num2str(num_run(2))];
  else  
    DEM_in    = double(Z_in);  
    file_out  = file_in;
  end
  log_file    = ['RESULTS/',file_out,'.log'];                               % Log file

  %... Get the coordinate of the area  
  [nrows,ncols] = size(DEM_in);      
  dx            = (bbox(2,1)-bbox(1,1))/Mx;
  dy            = (bbox(2,2)-bbox(1,2))/My;
  MAT.dx        = dx;
  MAT.dy        = dy;  
  
  xmin          = bbox(1,1);
  xmax          = xmin + ncols*dx;
  ymin          = bbox(1,2);
  ymax          = ymin + nrows*dy;
 
  %... Get current time for log file
  date_now = clock;                                                         
  if date_now(5) < 10
    date_sec = ['0',num2str(date_now(5))];
  else
    date_sec = num2str(date_now(5));
  end
  date_now = [num2str(date_now(2)),'/',num2str(date_now(3)),'/', ...
              num2str(date_now(1)),' at ',num2str(date_now(4)),':',date_sec];
  
  %... Open log file for printing. Also print to screen.
  fid = fopen(log_file,'at');
  disp(['Time: ',date_now]); 
  disp(['Study Area: ', file_in]); 
  
  fprintf(fid,['\n']);fprintf(fid,['\n']);fprintf(fid,['\n']);
  fprintf(fid,['Time: ',date_now,'\n']);
  fprintf(fid,['Study Area: ', file_in,'\n']);

  %... Convert negative (or missing data) into NaN
  nodata            = -998;    
  DEM_in(DEM_in<0)  = NaN;  
  
  %... DEM correction for filling data
  if (SWITCHES.correct)  
    disp(       ' DEM correction time   ');    
    fprintf(fid,' DEM correction time   \n');        
    tStart = tic;    
    if (interpolate_method == 1)
      % Convert missing data to one nodata value
      DEM_in(1,isnan(DEM_in(1,:)))      = min(min(DEM_in));
      DEM_in(end,isnan(DEM_in(end,:)))  = min(min(DEM_in));
      DEM_in(isnan(DEM_in(:,1)),1)      = min(min(DEM_in));
      DEM_in(isnan(DEM_in(:,end)),end)  = min(min(DEM_in));
      DEM_in(isnan(DEM_in))             = nodata;
      MAT.grid                          = DEM_in;  
  
      [TCA, Bound, Flow_dir] = mexD8(MAT.grid,MAT.dy/MAT.dx,flood,ndv);       % Get D8 flow direction
      [Conn_image, num_objects, Cen_multi, ind_mul, Cen_single, ind_sing]...  % Identify Multi-centers and Single-center for DEM correction
                    = D8_checkzero(Flow_dir, MAT.grid,mink);            
                
      % DEM corection using mex_DEM_correct to get nearest neighbor value 
      [MAT_cor]     = DEM_cor(MAT,Cen_multi,ind_mul,Cen_single,nodata,mink);  
    else
      % DEM correction using inpainting in image processing
      MAT.grid      = DEM_in;  
      MAT_cor.grid  = inpaint_nans(MAT.grid,3);      
      MAT_cor.dx    = dx;
      MAT_cor.dy    = dy;
    end
    
    MAT.grid      = MAT_cor.grid;                                           % Update DEM to corrected DEM
    MAT_ori.grid  = MAT_cor.grid;                                           % Back-up original DEM to corrected DEM        
    tEnd          = toc(tStart);    
    disp(['   Elapsed time is ', num2str(floor(tEnd/60)),' mins and ', ...
            num2str(rem(tEnd,60)),' secs']);
    fprintf(fid,'   Elapsed time is %d mins and %f secs\n', ...
            floor(tEnd/60),rem(tEnd,60));
  else
    disp(       ' No DEM correction   ');    
    fprintf(fid,' No DEM correction   \n');            
    MAT.grid = DEM_in;
  end
        
  %... Back-up original DEM and DEM at each level     
  MAT_ori = MAT;                                                      
  MAT_lev = MAT;    
  
  %... Write elevation data to to surfer format
  if (SWITCHES.surfer)
    namesurfer = ['RESULTS/',file_in,'.grd'];
    disp(['--------------------------------------------------']);
    disp(['      EXPORTING TO SURFER GRID FORMAT FILE        ']);
    disp(['--------------------------------------------------']);
    disp([' Bounding box: xmin = ',num2str(xmin),' (m)'       ]);
    disp(['               xmax = ',num2str(xmax),' (m)'       ]);
    disp(['               ymin = ',num2str(ymin),' (m)'       ]);
    disp(['               ymax = ',num2str(ymax),' (m)'       ]);
    disp(['                                                  ']);
    disp([' # of columns: ',num2str(ncols)                    ]);
    disp([' # of rows   : ',num2str(nrows)                    ]);
    disp(['                                                  ']);
    disp([' Grid size: dx = ',num2str(MAT.dx),' (m)'          ]);
    disp(['            dy = ',num2str(MAT.dy),' (m)'          ]);
    
    fprintf(fid,['--------------------------------------------------\n']);
    fprintf(fid,['      EXPORTING TO SURFER GRID FORMAT FILE        \n']);
    fprintf(fid,['--------------------------------------------------\n']);
    fprintf(fid,[' Bounding box: xmin = ',num2str(xmin),' (m)       \n']);
    fprintf(fid,['               xmax = ',num2str(xmax),' (m)       \n']);
    fprintf(fid,['               ymin = ',num2str(ymin),' (m)       \n']);
    fprintf(fid,['               ymax = ',num2str(ymax),' (m)       \n']);
    fprintf(fid,['                                                  \n']);
    fprintf(fid,[' # of columns: ',num2str(ncols),                 '\n']);
    fprintf(fid,[' # of rows   : ',num2str(nrows),                 '\n']);
    fprintf(fid,['                                                  \n']);
    fprintf(fid,[' Grid size: dx = ',num2str(MAT.dx),' (m)          \n']);
    fprintf(fid,['            dy = ',num2str(MAT.dy),' (m)          \n']);
    
    grd_write(MAT.grid,xmin,xmax,ymin,ymax,namesurfer);
  end
  
  
%% . . . MAIN PROGRAM CODE. . . . . . . . . . . . . . . . . . . . . . . . .
  %... Print statistics information  
  disp(['-------------------------------------------------'           ]);
  disp(['         DEPRESSION IDENTIFICATION MODEL         '           ]);
  disp(['-------------------------------------------------'           ]);
  disp([' Toporaphic data size: ', num2str(nrows),' x ',num2str(ncols)]);
  disp([' Grid size: dx = ',num2str(MAT.dx),' (m)'                    ]);
  disp(['            dy = ',num2str(MAT.dy),' (m)'                    ]);
  disp(['                                                 '           ]);
  disp([' Total # of searching process: ',num2str(ncols*nrows)        ]);
  disp(['                                                 '           ]);

  fprintf(fid,['-------------------------------------------------\n'  ]);
  fprintf(fid,['         DEPRESSION IDENTIFICATION MODEL         \n'  ]);
  fprintf(fid,['-------------------------------------------------\n'  ]);
  fprintf(fid,[' Toporaphic data size: ', num2str(nrows),' x ',num2str(ncols),'\n']);
  fprintf(fid,[' Grid size: dx = ',num2str(MAT.dx),' (m)\n'           ]);
  fprintf(fid,['            dy = ',num2str(MAT.dy),' (m)\n'           ]);
  fprintf(fid,['                                                 \n'  ]);
  fprintf(fid,[' Total # of searching process: ',num2str(ncols*nrows),'\n']);
  fprintf(fid,['                                                 \n'  ]);
  
  num_lev = input(' Enter the number of level: ');
  num_fig = 1;
          
  for lev = 1:num_lev
    disp( '                                                 ' );    
    disp( '                                                 ' );
    disp([' RUNNING LEVEL ',num2str(lev),'/',num2str(num_lev)]);
    disp( ' -------------------                             ' );    

    fprintf(fid, '                                                      \n' );    
    fprintf(fid, '                                                      \n' );
    fprintf(fid,[' SEARCHING LEVEL ',num2str(lev),'/',num2str(num_lev),'\n']);
    fprintf(fid, ' -------------------                                  \n' );  
  
  %. . .CENTERS & FLATS IDENTIFICATION. . . . . . . . . . . . . . . . . . .
    %..Identify centers (whose elevation is lower than 8 neighbors);
    %..If more than 1 cell have the same lowest elevation, they are a flat.
    disp(       ' CENTER SEARCH   '   );    
    fprintf(fid,' CENTER SEARCH   \n' );        
    
    tStart = tic;
    [TCA, Bound, Flow_dir] = mexD8(MAT.grid,MAT.dy/MAT.dx,flood,ndv);
    [Conn_image, num_objects, Cen_multi, ind_mul, Cen_single, ind_sing]...
          = D8_checkzero(Flow_dir, MAT.grid, mink);            
    tEnd  = toc(tStart);    
    
    disp(['   Elapsed time is ', num2str(floor(tEnd/60)),' mins and ', ...
          num2str(rem(tEnd,60)),' secs']);
    fprintf(fid,'   Elapsed time is %d mins and %f secs\n',...
          floor(tEnd/60),rem(tEnd,60));
 
 
  %. . .SEARCH & ADD: DEPRESSION & THRESHOLD CELLS. . . . . . . . . . . . .
    %..Load Center/Flat location list from D8_center function;
    %..A search process is activated from each center/flat;
    %..The process is stopped if threshold(s) are found or it reaches the 
    %     boundaries of the domain
    disp(       ['                      ']);    
    disp(       [' DEPRESSION SEARCH    ']);    
    fprintf(fid,['                      \n']);    
    fprintf(fid,[' DEPRESSION SEARCH    \n']);   
    
    tStart = tic;
    [Pud_pts,       Puddle,     ind_pud,    ind_pudcum,   ...
        Thres_pts,  Threshold,  ind_thres,  ind_threscum, ...
        volume,     area,       max_depth,  mean_depth]... 
      = pudsear(fid, MAT, MAT_ori, MAT_lev, Cen_multi, ind_mul, ...
                Cen_single, mink, max_size); 
    tEnd = toc(tStart);    
    
    disp(['   Elapsed time is ', num2str(floor(tEnd/60)),' mins and ', ...
          num2str(rem(tEnd,60)),' secs']);
    fprintf(fid,'   Elapsed time is %d mins and %f secs\n',...
          floor(tEnd/60),rem(tEnd,60));
    
  %. . .FILL UP DEPRESSIONS BY THRESHOLD ELEVATION FOR NEXT LEVEL SEARCH. .
    disp(       ['                     ']);    
    disp(       [' FILLING PROCESS     ']);    
    fprintf(fid,['                   \n']);        
    fprintf(fid,[' FILLING PROCESS   \n']);    
    
    tStart = tic;
    [DEM_pud] = mex_fill_pud(MAT.grid, double(Pud_pts)-1, ind_pudcum);
    MAT       = threshold_fill(MAT, DEM_pud,  Pud_pts,  ind_pudcum, ...
                  Thres_pts,ind_threscum,Threshold);    
    tEnd      = toc(tStart);    
    
    disp(['   Elapsed time is ', num2str(floor(tEnd/60)),' mins and ', ...
          num2str(rem(tEnd,60)),' secs']);
    fprintf(fid,'   Elapsed time is %d mins and %f secs\n',...
          floor(tEnd/60),rem(tEnd,60));
        
    MAT_lev   = MAT;
    
  %... Only show puddle whose area > 10 cells;  
    disp(       ['                        ']);    
    disp(       [' CONNECTING PROCESS     ']);  
    fprintf(fid,['                      \n']);            
    fprintf(fid,[' CONNECTING PROCESS   \n']);      
    tStart  = tic;
    PuddleN = Pud_connect(Puddle,10);
    tEnd    = toc(tStart);    
    
    disp(['   Elapsed time is ', num2str(floor(tEnd/60)),' mins and ', ...
            num2str(rem(tEnd,60)),' secs']);
    fprintf(fid,'   Elapsed time is %d mins and %f secs\n',...
            floor(tEnd/60),rem(tEnd,60));
    
  %% . . .PLOTTING. . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    x = 1:1:ncols;
    y = 1:1:nrows;
    if (SWITCHES.plots_on)
      lev_plotting;
    end

  %% . . .SAVING OUTPUTS. . . . . . . . . . . . . . . . . . . . . . . . . .
    if (SWITCHES.save_on)
      filename = ['./RESULTS/',file_out,'_LEV_',num2str(lev),'.mat'];      
      save(filename,'Puddle',     'Pud_pts',    'ind_pud',    'PuddleN',    ...
                    'Threshold',  'Thres_pts',  'ind_thres',  'Flow_dir',   ...
                    'Cen_multi',  'Cen_single', 'ind_sing',   'ind_mul',    ...
                    'area',       'volume',     'max_depth',  'mean_depth', ... 
                    'xmax',       'xmin',       'ymax',       'ymin',       ...
                    'nrows',      'ncols');
      filetiff = ['./RESULTS/',file_out,'_LEV_',num2str(lev),'.tif'];       
      geotiffwrite(filetiff,single(Puddle),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    end
  end
    
%% . . .FINAL DISPLAY. . . . . . . . . . . . . . . . . . . . . . . . . . . 
  disp('                                                     ');
  disp(' COMPLETED!!!                                        ');
  disp('                                                     ');
  disp('                                                     ');
  disp('                                                     ');
  fprintf(fid,'\n'                                            );
  fprintf(fid,'. . . . . . . . . COMPLETED!!! . . . . . . . .');
  fprintf(fid,'\n'                                            );
  fprintf(fid,'\n'                                            );
  fprintf(fid,'\n'                                            );

%% . . .CLOSE WRITING LOG FILE. . . . . . . . . . . . . . . . . . . . . . . 
  fclose(fid);                                                 
