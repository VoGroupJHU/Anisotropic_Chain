function [wavefunction_full,H_full,indx_zr] = compute_wavefunction_core(pts_lattice,q_lattice,type_lattice,coord_lattice,sigma_core,Fk_lattice,rscale_lattice,box_shift,grid_step)
% Inputs
% pts_lattice:      lattice position
% q_lattice:        lattice orientation
% type_lattice:     lattice polyhedra type
% coord_lattice:    vertices shape
% Fk_lattice:       parameterized shape kernel base shape
% Fk_patch_lattice: parameterized shape kernel patchy face
% patch_core_ratio: strength of patch relative to core
% Fk_bond_lattice:  parameterized shape kernel bond
% bond_core_ratio:  strength of bond relative to core
% rscale_lattice:   orbital scaling ratio
% box_shift:        box expansion buffer factor
% grid_step:        step size of meshgrid

% Outputs
% wavefunction_full:    evaluated wavefunction
% H_full:               evaluated Hamiltonian
% indx_zr:              index of points inside polyhedra 

% Generate full grid
[xspace,yspace,zspace,pts_grid,xx] = generate_meshgrid(pts_lattice,box_shift,grid_step);

% Determine interior point
indx_interior = zeros(length(pts_grid(:,1)),1);
for i = 1:length(pts_lattice(:,1))
    % Define polyhedra coordinates
    pts_coord = coord_lattice{type_lattice(i)};
    % Rotate and shift
    tmp = pts_coord*rot_quat(q_lattice(i,:));
    tmp = bsxfun(@plus,tmp,pts_lattice(i,:));
    % Create alphaShape
    shape_tmp = alphaShape(tmp,1E10);
    % Determine interior points
    indx_tmp = inShape(shape_tmp,pts_grid);
    % Update master list
    indx_interior = indx_interior + indx_tmp;
end
indx_zr = find(indx_interior ~= 0);

% Calculate wavefunction
n_lv = 1.0;
H_full = zeros(size(xx));
wave_full = zeros(size(xx));
wavefunction_full = cell(1,length(pts_lattice(:,1)));
for i = 1:length(pts_lattice(:,1))
    % Grab lattice point
    pts_tmp = pts_lattice(i,:);
    % Scaling distance
    rscale = rscale_lattice(type_lattice(i));
    % Rotate to body frame and determine distances
    dr_tmp = bsxfun(@minus,pts_grid,pts_tmp);
    q_tmp = q_lattice(i,:);
    q_tmp(2:end) = -q_tmp(2:end);
    dr_tmp = dr_tmp*rot_quat(q_tmp);
    % Convert to angle
    [theta,phi,r_tmp] = xyz_to_kernel(dr_tmp);
    % Define and compute kernel
    Fk = Fk_lattice{type_lattice(i)};
    kernel_tmp = Fk(theta,phi);
    % Calculate wavefunction
    Kpp = 0.5*sqrt(4*4*3*kernel_tmp + 1) - 0.5;
    % Base shape wavefunction
    V_tmp = -2*kernel_tmp.^(2).*(r_tmp*rscale)/n_lv; 
    wave_tmp = bsxfun(@power,r_tmp,Kpp).*exp(V_tmp);
    wave_tmp(indx_zr) = 0;
    % Integrate and normalize
    wave_int = reshape(wave_tmp,size(xx));
    I = trapz(zspace/rscale,trapz(xspace/rscale,trapz(yspace/rscale,wave_int.^2)));
    wave_int = wave_int/sqrt(I);
    % Add
    wave_full = wave_full + wave_int;
    H_full = H_full - reshape(kernel_tmp,size(xx))./reshape(r_tmp,size(xx));
    H_full(indx_zr) = 0;
    % Store individual wavefunctions
    wavefunction_full{i} = wave_int;
end
% Core-core repulsion
for i = 1:length(pts_lattice(:,1))
    % Grab orientation
    q1 = q_lattice(i,:);
    for j = i:length(pts_lattice(:,1))
        % Grab orientation
        q2 = q_lattice(i,:);
        % Grab non-identical particle
        if i ~= j
            % Normal vector
            dr = pts_lattice(i,:) - pts_lattice(j,:);
            k12 = 1;
            H_full = H_full + (k12*sigma_core./norm(dr)).^12 - 1*(k12*sigma_core./norm(dr)).^6;
        end
    end
end

end