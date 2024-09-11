function rpeak_scale = calc_rscale(Fk)
% Define ground state constants
p_power = 0;
n_lv = 1.0;

% Generate meshgrid
Lx = 15;
Ly = 15;
Lz = 15;
dstep = 0.1;
xspace = -Lx:dstep:Lx;
yspace = -Ly:dstep:Ly;
zspace = -Lz:dstep:Lz;
[xx,yy,zz] = meshgrid(xspace,yspace,zspace);

% Convert to grid
pts_grid = [xx(:) yy(:) zz(:)];

% Convert to angle
[theta,phi,r_tmp] = xyz_to_kernel(pts_grid);
kernel_tmp = Fk(theta,phi);
% Calculate wavefunction
Kpp = 0.5*sqrt(4*4*3*kernel_tmp + 1) - 0.5;
wave_tmp = bsxfun(@power,r_tmp,Kpp).*exp(-2*kernel_tmp.^(2-p_power).*r_tmp/n_lv);
wave_int = reshape(wave_tmp,size(xx));
I = trapz(zspace,trapz(yspace,trapz(xspace,wave_int.^2)));
wave_int = wave_int/sqrt(I);
% Determine scaling distance for peak
[~,indx_scale] = max(wave_int(:));
rpeak_scale = r_tmp(indx_scale);
end