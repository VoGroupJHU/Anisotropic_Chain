function [wavefunction_full,H_full,indx_zr] = compute_wavefunction_bond(pts_lattice,q_lattice,type_lattice,coord_lattice,sigma_core,Fk_lattice,rscale_lattice,pts_bond_lattice,sigma_bond_lattice,bond_list,box_shift,grid_step)
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
    %%% Loop over bonds %%%
    pts_bond = pts_bond_lattice{type_lattice(i)};
    sigma_bond = sigma_bond_lattice{type_lattice(i)};
    for ii = 1:length(pts_bond(:,1))
        % Rotate and shift bond location
        tmp_bond = pts_bond(ii,:)*rot_quat(q_lattice(i,:));
        tmp_bond = bsxfun(@plus,tmp_bond,pts_lattice(i,:));
        % Compute distance
        dr_bond = bsxfun(@minus,pts_grid,tmp_bond);
        % Check distance
        dr_bond = sqrt(sum(dr_bond.^2,2));
        indx_bond = dr_bond < (sigma_bond(ii)/2);
        % Update bond
        indx_interior = indx_interior + indx_bond;
    end
    %%%%%%%%%%%%%%%%%%%%%%%    
end
indx_zr = find(indx_interior ~= 0);

% Calculate wavefunction
n_lv = 1.0;
H_full = zeros(size(xx));
wave_full = zeros(size(xx));
wavefunction_full = cell(1,length(pts_lattice(:,1)));
for i = 1:length(pts_lattice(:,1))
    % disp([i length(pts_lattice(:,1))])
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
    % Define scaling
    wave_base_max = max(abs(wave_tmp(:)));
    %%% Loop over bonds %%%    
    wave_bond_total = 0*wave_tmp;
    pts_bond = pts_bond_lattice{type_lattice(i)};
    sigma_bond = sigma_bond_lattice{type_lattice(i)};
    for ii = 1:length(pts_bond(:,1))
        % Grab lattice point
        pts_bond_tmp = pts_lattice(i,:);
        % Shift to bond center
        dr_bond = bsxfun(@minus,pts_grid,pts_bond_tmp);        
        % Rotate to body frame and determine distances
        q_tmp = q_lattice(i,:);
        q_tmp(2:end) = -q_tmp(2:end);
        dr_bond = dr_bond*rot_quat(q_tmp);
        % Shift
        dr_bond = bsxfun(@minus,dr_bond,pts_bond(ii,:));
        % Compute distance
        r_bond = sqrt(sum(dr_bond.^2,2));
        % Kernel
        kernel_bond = 1;
        % Diameter scaling
        sigma_scale = (0.5*sigma_bond(ii)).^(-1);
        % Calculate wavefunction
        Kpp = 0.5*sqrt(4*4*3*kernel_bond + 1) - 0.5;
        % Bond potential term
        bond_tmp = sigma_bond(ii);
        V_bond = -2*kernel_bond.^(2).*(r_bond*rscale/(bond_tmp))/n_lv;
        wave_bond = bsxfun(@power,r_tmp,Kpp).*exp(V_bond);
        % Scale
        wave_bond = wave_bond./max(abs(wave_bond(:)))*wave_base_max;
        % Update
        wave_bond_total = wave_bond_total + wave_bond;
    end
    % Zeros
    wave_bond_total(indx_zr) = 0;
    wave_tmp(indx_zr) = 0;
    % Combine
    wave_tmp = 0*wave_tmp + wave_bond_total;
    %%%%%%%%%%%%%%%%%%%%%%%    
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
% Bond potential
for i = 1:length(pts_lattice(:,1))
    pts_bond_i = pts_bond_lattice{type_lattice(i)};
    for j = i:length(pts_lattice(:,1))
        pts_bond_j = pts_bond_lattice{type_lattice(j)};
        % Check bond list
        flag_bond = 0;
        for k = 1:length(bond_list(:,1))
            bond_tmp = [i j];
            bond_tmp = unique(bond_tmp);
            bond_check = isequal(bond_tmp,unique(bond_list(k,1:2)));
            if bond_check == 1
                flag_bond = 1;
                indx_bond_list = k;
                break
            end
        end
        % Compute bond potential
        if flag_bond == 1
            % Loop through all bonding locations            
            for ii = 1:length(pts_bond_i(:,1))
                % Place first bond
                tmp_bond1 = pts_bond_i(ii,:)*rot_quat(q_lattice(i,:));
                tmp_bond1 = bsxfun(@plus,tmp_bond1,pts_lattice(i,:));
                for jj = 1:length(pts_bond_j(:,1))
                    % Rotate and shift bond location                    
                    % Place second bond
                    tmp_bond2 = pts_bond_j(jj,:)*rot_quat(q_lattice(j,:));
                    tmp_bond2 = bsxfun(@plus,tmp_bond2,pts_lattice(j,:));
                    % Distance
                    dr_bond = norm(tmp_bond2-tmp_bond1);                                    
                    % Add bond energy
                    k_bond = bond_list(k,5);
                    dr_bond0 = bond_list(k,6);
                    % Harmonic
                    if ( (bond_list(indx_bond_list,3) == ii) && (bond_list(indx_bond_list,4) == jj) )
                        H_full = H_full + 0.5*k_bond*(dr_bond - dr_bond0).^2;
                    else
                        H_full = H_full + 0;
                    end
                end
            end
        end
    end
end
end