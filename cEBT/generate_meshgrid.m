function [xspace,yspace,zspace,pts_grid,xx,yy,zz] = generate_meshgrid(pts_lattice,box_shift,grid_step)

% Generate full grid
Lx = max(abs(pts_lattice(:,1))) + box_shift;
Ly = max(abs(pts_lattice(:,2))) + box_shift;
Lz = max(abs(pts_lattice(:,3))) + box_shift;

% Define box size
LT = max([Lx Ly Lz]);
Lx = LT;
Ly = LT;
Lz = LT;

% Generate meshgrid
dstep = grid_step;
xspace = -Lx:dstep:Lx;
yspace = -Ly:dstep:Ly;
zspace = -Lz:dstep:Lz;
[xx,yy,zz] = meshgrid(xspace,yspace,zspace);

% Convert to grid
pts_grid = [xx(:) yy(:) zz(:)];