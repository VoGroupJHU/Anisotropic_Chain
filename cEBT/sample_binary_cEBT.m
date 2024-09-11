clear; close all; clc;

warning off;
addpath('C:\Users\lanth\Documents\codes\polyhedron_function')

%%%%%%%%%%%%%%%%%%%%%
%%% Define colors %%%
%%%%%%%%%%%%%%%%%%%%%

purple = [107.0/255, 62.0/255, 154.0/255];
blue = [0, 64.0/255, 128.0/255];
blue2 = [85 170 255]/255;
red = [180.0/255, 7.0/255, 5.0/255];
red2 = [206 0 0 ]/206;
green = [0.0/255, 112/255, 0.0/255];
gold = [197 190 25]/255;
gray = [137 137 137]/255;
gray2 = [247 247 247]/255;
orbital_color = [171 202 247]/255;

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
%%% Define colormap %%%
%%%%%%%%%%%%%%%%%%%%%%%

% Define colormap
% Blue
color0 = [0 25 51];
color1 = [0 76 153];
color2 = [0 102 204];
color3 = [0 128 255];
color4 = [51 153 255];
color5 = [102 179 255];
color6 = [153 204 255];
color7 = [192 192 192];
color8 = [224 224 224];
color_mtx = [color7; color6; color5; color4; color3; color2; color1]/255;
% Generate grid
n_color = 1000;
n_color_type = length(color_mtx(:,1));
n_color_grid = ceil(n_color/n_color_type);
mymap = [];
for i = 1:n_color_type-1
    tmp_ith = color_mtx(i,:);
    tmp_jth = color_mtx(i+1,:);
    tmp_mtx = [linspace(tmp_ith(1),tmp_jth(1),n_color_grid)', linspace(tmp_ith(2),tmp_jth(2),n_color_grid)', linspace(tmp_ith(3),tmp_jth(3),n_color_grid)'];
    mymap = [mymap; tmp_mtx];
end

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate and define shape/kernel %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load and define shape
shape_name = 'Cube';
pts_coord = textread(strcat(shape_name,'.txt'));
pts_coord = bsxfun(@minus,pts_coord,mean(pts_coord));

% Generate polyhedra
grid_point = 20;
pts_poly = gen_poly(pts_coord,grid_point);
pts_poly = roundn(pts_poly,-3);
pts_poly = unique(pts_poly,'rows');

% Insphere and circumsphere
dr = sqrt(sum(pts_poly.^2,2));
R_insph = min(dr(:));
R_circum = max(dr(:));

% Particle size
sigma_core = 2*R_insph;
sigma_core_out = 2*R_circum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define faces, edges, corners %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Faces %%%
% Define faces
if strcmp(shape_name,'Cube') == 1
    faces = {[1 2 6 5],[5 6 8 7],[3 7 8 4],[3 4 2 1],[5 1 3 7],[2 6 8 4]};
elseif strcmp(shape_name,'TruncatedCube') == 1
    faces = {[21, 2, 20], [19, 1, 22], [3, 17, 24], [23, 18, 4], ...
        [5, 11, 14], [16, 9, 7], [13, 12, 6], [8, 10, 15], ...
        [18, 23, 24, 17, 19, 22, 21, 20], ...
        [7, 9, 11, 5, 1, 19, 17, 3], [2, 21, 22, 1, 5, 14, 13, 6], ...
        [8, 15, 16, 7, 3, 24, 23, 4], [12, 13, 14, 11, 9, 16, 15, 10], ...
        [4, 18, 20, 2, 6, 12, 10, 8]};
elseif strcmp(shape_name,'TruncatedTetrahedron') == 1
    faces = {[5 7 6], [6 7 8 12 2 4], [6 4 10 9 3 5], [1 9 3], [12 8 11],...
        [7 5 3 1 11 8], [10 4 2]};
elseif strcmp(shape_name,'TruncatedOctahedron') == 1
    faces = {[22 24 21 18],[6 12 5 2],[15 17 11 9],[14 16 10 8],[23 20 13 19],...
        [7 4 1 3],[14 19 13 7 3 8],[18 21 16 10 5 12],[17 22 18 12 6 11],...
        [9 11 6 2 1 4],[4 7 13 20 15 9],[15 20 23 24 22 17],[19 14 16 21 24 23],...
        [1 2 5 10 8 3]};
elseif strcmp(shape_name,'RhombicDodecahedron') == 1
    faces = {[7, 5, 1, 2], [12, 11, 5, 7], [2, 4, 9, 7], [7, 9, 14, 12], ...
        [2, 1, 3, 4], [6, 1, 5, 11], [3, 1, 6, 8], [14, 9, 4, 10], ...
        [8, 10, 4, 3], [13, 11, 12, 14], [8, 6, 11, 13], [13, 14, 10, 8]};
elseif strcmp(shape_name,'Dodecahedron') == 1
    faces = {[8 4 15 1 16],[7 16 1 14 3],[1 15 10 9 14],[7 11 12 8 16],...
        [2 6 12 11 5],[5 11 7 3 19],[19 3 14 9 17],[9 10 18 13 17],...
        [6 12 8 4 20],[20 4 15 10 18],[2 13 18 20 6],[2 5 19 17 13]};
elseif strcmp(shape_name,'Prism10') == 1
    faces = {[1 7 11 15 19 3 17 13 9 5],[2 8 12 16 20 4 18 14 10 6],...
        [6 10 9 5],[10 14 13 9],[14 18 17 13],[18 4 3 17],[4 20 19 3],...
        [20 16 15 19],[16 12 11 15],[12 8 7 11],[8 2 1 7],[2 6 5 1]};     
else
    faces_tmp = boundary(pts_coord,0);
    faces = {};
    for i = 1:length(faces_tmp(:,1))
        faces{i} = faces_tmp(i,:);
    end
end
%%%%%%%%%%%%%

%%% Edges %%%
edges_tmp = [];
% Pairwise distance
pdist_edge = pdist2(pts_coord,pts_coord);
% Round
pdist_edge = roundn(pdist_edge,-3);
% Remove zeros
indx_zr = find(pdist_edge == 0);
pdist_edge(indx_zr) = 1E10;
% Loop and detect distance
for i = 1:length(pdist_edge(:,1))
    % Row
    pdist_tmp = pdist_edge(i,:);
    % Shortest distance
    pdist_tmp_min = min(pdist_tmp);
    % Index of points with shortest distance
    indx_tmp = find(pdist_tmp == pdist_tmp_min);
    % Add to edge
    for ii = 1:length(indx_tmp)
        edge_add = [i indx_tmp(ii)];
        edges_tmp = [edges_tmp; edge_add];
    end
end
% Get unique edges
edges_tmp = sort(edges_tmp,2);
edges_tmp = unique(edges_tmp,'rows');
% Define final
edges = cell(1,length(edges_tmp(:,1)));
for i = 1:length(edges_tmp(:,1))
    edges{i} = edges_tmp(i,:);
end
%%%%%%%%%%%%%

%%% Corners %%%
corners = pts_coord;
%%%%%%%%%%%%%%%

%%% Define bond %%%
pts_bond1 = pts_coord(3,:);
pts_bond = [pts_bond1];
%%%%%%%%%%%%%%%%%%%

%%% Align bond along z-axis %%%
nv_rot = pts_bond/norm(pts_bond);
pts_coord = pts_coord*rot_az(nv_rot,1);
pts_bond = pts_bond*rot_az(nv_rot,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define bond space %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Define bond
bond_size_array = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loop through bond size array %%%
for bbb = 1:length(bond_size_array)

    %%% Define full bond set %%%
    bond_size = (bond_size_array(bbb)*R_insph)*ones(length(pts_bond1(:,1)),1);
    sigma_bond = bond_size*ones(length(pts_bond(:,1)),1);
    % Bond list
    % Order: ith index; jth index; bond_location_indx_ith; bond_location_indx_jth; bond_spring; bond_length
    k_bond = 30;
    bond_length = 1.5*sigma_bond(1);
    bond_list = [1 2 1 1 k_bond bond_length];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Define orientational space %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Angle step
    theta_step = 15;
    phi_step = theta_step;

    % Define range theta
    theta_start = -90;
    theta_end = 90;

    % Define range phi
    phi_start = -90;
    phi_end = 90;

    % Define space for theta
    theta1 = theta_start:2*theta_step:theta_end;
    if roundn(theta1(end),-2) ~= roundn(theta_end,-2)
        theta1(end+1) = theta_end;
    end
    % Add zero if needed
    if isempty(intersect(roundn(theta1,-2),0))
        theta1(end+1) = 0;
    end
    % Sort
    theta1 = sort(theta1);

    % Define space for phi
    phi1 = phi_start:phi_step:phi_end;
    if roundn(phi1(end),-2) ~= roundn(phi_end,-2)
        phi1(end+1) = phi_end;
    end
    % Add zero if needed
    if isempty(intersect(roundn(phi1,-2),0))
        phi1(end+1) = 0;
    end

    % Second particle
    theta1 = 0;
    phi1 = 0;
    theta2 = theta1 - 180;
    phi2 = phi1 - 0*180;
    % phi2 = phi1 - 180;

    % Define phase space
    phase_space = zeros(numel(theta1)^2*numel(phi1)^2,4);
    counter = 1;
    for n = 1:length(theta1)
        for nn = 1:length(phi1)
            for m = 1:length(theta2)
                for mm = 1:length(phi2)
                    % Define angle
                    tmp = [theta1(n) phi1(nn) theta2(m) phi2(mm)];
                    % Convert to radians
                    tmp = tmp*pi/180;
                    phase_space(counter,:) = tmp;
                    counter = counter + 1;
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Generate kernel %%%
    %%%%%%%%%%%%%%%%%%%%%%%

    %%% Core kernel %%%
    dgrid = 20;
    Fk_core = gen_kernel_evaluator(pts_coord,dgrid);
    %%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Determine distance scaling for single particle %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p_power = 0.0;
    n_lv = 1.0;
    rscale = calc_rscale(Fk_core);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Loop over phase space %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Meshgrid parameters
    box_shift = 3;
    grid_step = 0.1;

    % Initialize
    E_phase_space = zeros(length(phase_space(:,1)),1);
    E_full = cell(length(phase_space(:,1)),1);
    r_full = cell(length(phase_space(:,1)),1);
    pts_lattice_full = cell(length(phase_space(:,1)),1);

    % Define parallel
    % parpool(4)

    for pp = 1:length(phase_space(:,1))

        disp([pp length(phase_space(:,1)) bbb length(bond_size_array)])


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Compute wavefunctions and compute energy %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for nn = 1:2

            %%% Reset separation grid %%%
            rspace = linspace(1,3,100);
            rspace = 0.75:0.05:3;
            % rspace = 0.75:0.01:3;
            rspace = 0.65:0.05:3;
            rspace_plot = rspace;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%% Identify minimum from second loop %%%
            if nn == 2
                % Get minimum energy
                [Emin,indx_min] = min(Espace);
                % Store
                E_phase_space(pp) = Emin;
                % Store value for plotting
                E_full{pp} = Espace;
                r_full{pp} = rspace;
                rspace = rspace(indx_min);
                rspace = 1.4;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%% Energy calculation %%%
            % Initialize
            Espace = zeros(size(rspace));
            % Loop through separation grid
            for rr = 1:length(rspace)

                %%% Define lattice %%%
                % Place particles
                pts1 = [0 0 0];
                pts2 = [0 0 rspace(rr)*(2*R_circum)];
                % Shift and center
                if nn == 1
                    pts_lattice = [pts1; pts2];
                else
                    pts_lattice = pts1;
                end
                pts_lattice = [pts1; pts2];
                pts_lattice = bsxfun(@minus,pts_lattice,mean(pts_lattice));
                % Lattice type
                type_lattice = [1 1];
                %%% Orientation %%%
                % First particle
                q1_theta = [cos(phase_space(pp,1)/2) 0 0 sin(phase_space(pp,1)/2)];
                q1_phi = [cos(phase_space(pp,2)/2) 0 sin(phase_space(pp,2)/2) 0 ];
                rtx1 = rot_quat(q1_theta)*rot_quat(q1_phi);
                q1 = rotm2quat(rtx1);
                % Second particle
                q2_theta = [cos(phase_space(pp,3)/2) 0 0 sin(phase_space(pp,3)/2)];
                q2_phi = [cos(phase_space(pp,4)/2) 0 sin(phase_space(pp,4)/2) 0 ];
                rtx2 = rot_quat(q2_theta)*rot_quat(q2_phi);
                q2 = rotm2quat(rtx2);
                % Lattice
                q_lattice = [q1; q2];
                %%%%%%%%%%%%%%%%%%%
                % Vertices of coordinates for lattice
                coord_lattice = {pts_coord};
                % Scaling distance
                rscale_lattice = [rscale];
                % Core kernel
                Fk_core_lattice = {Fk_core};
                % Bond location
                pts_bond_lattice = {pts_bond};
                % Size of patch
                sigma_bond_lattice = {sigma_bond};
                % Energy ratios
                bond_energy = 1;
                %%%%%%%%%%%%%%%%%%%%%%

                %%% Generate meshgrid %%%
                if nn == 2
                    grid_step = 0.1;
                else
                    grid_step = 1/3;
                end
                [xspace,yspace,zspace,pts_grid,xx,yy,zz] = generate_meshgrid(pts_lattice,box_shift,grid_step);
                %%%%%%%%%%%%%%%%%%%%%%%%%

                %%% Individual wavefunction %%%
                % sigma_core = sigma_core_out;
                [wavefunction_core,H_core,indx_zr] = compute_wavefunction_core(pts_lattice,q_lattice,type_lattice,coord_lattice,sigma_core_out,Fk_core_lattice,rscale_lattice,box_shift,grid_step);
                [wavefunction_bond,H_bond] = compute_wavefunction_bond(pts_lattice,q_lattice,type_lattice,coord_lattice,sigma_core_out,Fk_core_lattice,rscale_lattice,pts_bond_lattice,sigma_bond_lattice,bond_list,box_shift,grid_step);
                % Create composite wavefunction
                wavefunction_full = cell(1,length(pts_lattice(:,1)));
                for i = 1:length(pts_lattice(:,1))
                    % Normalize
                    wave_core_add = wavefunction_core{i}/max(abs(wavefunction_core{i}(:)));
                    wave_bond_add = wavefunction_bond{i}/max(abs(wavefunction_bond{i}(:)));
                    % Combine wavefunctions
                    wave_full = wave_core_add + bond_energy*wave_bond_add;
                    wavefunction_full{i} = wave_full;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%% Solve for energy %%%
                % Core energy
                E_core = compute_energy(pts_lattice,wavefunction_core,H_core,xspace,yspace,zspace);
                % Bond energy
                E_bond = compute_energy(pts_lattice,wavefunction_bond,H_bond,xspace,yspace,zspace);
                % Combined
                E_tmp = (E_core + bond_energy*E_bond)/(1+bond_energy);
                % Store
                Espace(rr) = E_tmp;
                %%%%%%%%%%%%%%%%%%%%%%%%


            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%%% Plotting %%%
%%%%%%%%%%%%%%%%

rscale_shape = 1;

% Combined
figure(3)
wave_core = zeros(size(xx));
for i = 1:length(wavefunction_full)
    wave_core = wave_core + wavefunction_full{i};
end
I = trapz(zspace,trapz(yspace,trapz(xspace,wave_core.^2)));
wave_core = wave_core/sqrt(I);
% Identify isosurface
isovalue = mean(wave_core(:)) + 4.5*std(wave_core(:));
hold on;
fv = isosurface(xx,yy,zz,wave_core,isovalue);
p = patch(fv,'FaceColor',orbital_color,'EdgeColor','none','FaceAlpha',0.5);
isonormals(xx,yy,zz,wave_core,p)
lighting gouraud
%%% Orientation %%%
% First particle
q1_theta = [cos(phase_space(1,1)/2) 0 0 sin(phase_space(1,1)/2)];
q1_phi = [cos(phase_space(1,2)/2) 0 sin(phase_space(1,2)/2) 0 ];
rtx1 = rot_quat(q1_theta)*rot_quat(q1_phi);
q1 = rotm2quat(rtx1);
% Second particle
q2_theta = [cos(phase_space(1,3)/2) 0 0 sin(phase_space(1,3)/2)];
q2_phi = [cos(phase_space(1,4)/2) 0 sin(phase_space(1,4)/2) 0 ];
rtx2 = rot_quat(q2_theta)*rot_quat(q2_phi);
q2 = rotm2quat(rtx2);
% Lattice
q_lattice = [q1; q2];
%%%%%%%%%%%%%%%%%%%
for i = 1:length(pts_lattice(:,1))
    tmp = bsxfun(@plus,(1.1*rscale_shape)*pts_coord*rot_quat(q_lattice(i,:)),pts_lattice(i,:));
    tmp_shp = alphaShape(tmp,10);
    plot_core = plot(tmp_shp,'FaceColor',gray,'Edgecolor','none');
    if length(faces) > 0
        for j = 1:length(faces)
            k = faces{j};
            k = [k k(1)];
            plot3(tmp(k,1),tmp(k,2),tmp(k,3),'k','Linewidth',3.5)
        end
    end
    set(plot_core,'facelighting','flat')
end
view_angle = [21 25];
camlight
axis equal
camlight(view_angle(1)-30,view_angle(2)+60)
view(view_angle)
xlabel('x')
ylabel('y')
zlabel('z')
rotate3d('on')
axis off

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%