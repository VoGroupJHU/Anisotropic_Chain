function E_compute = compute_energy(pts_lattice,wavefunction_full,H_full,xspace,yspace,zspace)

%%% Solve for energy %%%
% Solve for energy
mtx_S = zeros(length(pts_lattice(:,1)));
mtx_H = zeros(length(pts_lattice(:,1)));
for i = 1:length(pts_lattice(:,1))
    disp([i length(pts_lattice(:,1))])
    for j = 1:length(pts_lattice(:,1))        
        wave_tmp = wavefunction_full{i}.*wavefunction_full{j};
        mtx_S(i,j) = trapz(zspace,trapz(xspace,trapz(yspace,wave_tmp)));
        wave_tmp = wavefunction_full{i}.*H_full.*wavefunction_full{j};
        mtx_H(i,j) = trapz(zspace,trapz(xspace,trapz(yspace,wave_tmp)));
    end
end
% Construct matrix
mtx_eff = inv(mtx_S)*mtx_H;
[V,D] = eig(mtx_eff);
E_compute = min(diag(D));
%%%%%%%%%%%%%%%%%%%%%%%%