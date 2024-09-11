clear; close all; clc;

warning off;

% Paper: "Elucidating the Interplay Between Entropy-Driven and Patch-Mediated 
% Bonding in Directing Nanoscale Assemblies"
% Sample Code for chain initialization
% Authors: Kireeti Akkunuri, Xiangyu Zhang, Thi Vo
% Dated: September 2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define system parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARAMETER -- Chain length
N = 40;


% PARAMETER -- monomer shape
shape_name = 'Tetrahedron';


%%%%%%%%%%%%%%%%%%%%
%%% Scale system %%%
%%%%%%%%%%%%%%%%%%%%

% Load shape from file with vertices
%% Vertices in shape file correspond to shape with insphere diameter of 1
pts_verts = textread(strcat(shape_name,'.txt'));
pts_verts = bsxfun(@minus,pts_verts,mean(pts_verts));

% Define faces
if strcmp(shape_name,'Cube') == 1
    faces = {[1 2 6 5],[5 6 8 7],[3 7 8 4],[3 4 2 1],[5 1 3 7],[2 6 8 4]};
else    
    faces_tmp = boundary(pts_verts,0);
    faces = {};
    for i = 1:length(faces_tmp(:,1))
        faces{i} = faces_tmp(i,:);
    end
end

% Determine kernel
grid_poly = 100;
pts_poly = gen_poly(pts_verts,grid_poly);       % use the prototype vertices (centered about 000)
dr_poly = sqrt(sum(pts_poly.^2,2));
% Insphere - for scaling
r_insph = min(dr_poly(:));  % smallest inspher distance from any of the grid points

% Circumsphere
r_circum = max(dr_poly(:));

% Average sigma
r_avg = 0.5*(r_insph+r_circum);
% r_avg = 0.75*r_insph + 0.25*r_circum;

% Define ghost (=patch) size
sigma_ghost = 10.0;
v_ghost = (4/3)*pi*(sigma_ghost/2)^3;

% Define ratio of ghost to monomer insphere
R_ghost_monomer = 2;


% Vertics scale factor
verts_scale = R_ghost_monomer/r_insph*sigma_ghost;
verts_scale = roundn(verts_scale,-6);

% Scale
pts_verts_scaled = verts_scale*pts_verts;

% For tetrahedron
sigma_ghost = 10.0*(3/1.732)
% For cube and octahedron, sigma_ghost = 10 


% Scale radii
r_insph = verts_scale*r_insph;
r_circum = verts_scale*r_circum;
r_avg = verts_scale*r_avg;

% Compute volume
[~,v_shape] = boundary(pts_verts_scaled,0);

% Compute moment of inertia
I_tensor = calc_inertial_tensor(pts_verts_scaled);

%%% Define edge %%%
edges = [];
% Pairwise distance between rows/points
pdist_edge = pdist2(pts_verts_scaled,pts_verts_scaled);
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
        edges = [edges; edge_add];
    end
end
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Place ghost particle %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%e

% Scale parameter
ghost_scale_face = 1.0;
ghost_scale_edge = 1.00;
ghost_scale_corner = 1.05;
ghost_scale_corner = 1.0;


%%% Face ghost %%%
face_center = zeros(length(faces),3);
for i = 1:length(faces)
    ktmp = faces{i};
    pts_face = pts_verts_scaled(ktmp,:);
    face_center(i,:) = mean(pts_face);
end
% First ghost
pts_ghost1 = face_center(1,:);
% Determine position of ghost
dr_ghost1 = bsxfun(@minus,face_center,pts_ghost1);
dr_ghost1 = sqrt(sum(dr_ghost1.^2,2));
[~,indx_ghost2] = max(dr_ghost1);
% Second ghost
pts_ghost2 = face_center(indx_ghost2,:);
% Place ghost
pts_face_ghost = [pts_ghost1; pts_ghost2];
% % Define scaling
% for i = 1:length(pts_face_ghost(:,1))
%     pts_tmp = pts_face_ghost(i,:);
%     nv_tmp = pts_tmp/norm(pts_tmp);
%     pts_tmp = pts_tmp + nv_tmp*(0.15*(r_insph+sigma_ghost/2));
%     pts_face_ghost(i,:) = pts_tmp;
% end
% Scale
pts_face_ghost = ghost_scale_face*pts_face_ghost;
%%%%%%%%%%%%%%%%%%

%%% Corner ghost %%%
% First ghost
pts_ghost1 = pts_verts_scaled(1,:);
% Determine position of ghost
dr_ghost1 = bsxfun(@minus,pts_verts_scaled,pts_ghost1);
dr_ghost1 = sqrt(sum(dr_ghost1.^2,2));
[~,indx_ghost2] = max(dr_ghost1);
% Second ghost
pts_ghost2 = pts_verts_scaled(indx_ghost2,:);
% Place ghost
pts_corner_ghost = [pts_ghost1; pts_ghost2];
% % Define scaling
% for i = 1:length(pts_corner_ghost(:,1))
%     pts_tmp = pts_corner_ghost(i,:);
%     nv_tmp = pts_tmp/norm(pts_tmp);
%     pts_tmp = pts_tmp + nv_tmp*(0.15*(r_insph+sigma_ghost/2));
%     pts_corner_ghost(i,:) = pts_tmp;
% end
% Scale
pts_corner_ghost = ghost_scale_corner*pts_corner_ghost;
%%%%%%%%%%%%%%%%%%%%

%%% Edge ghost %%%
% Edge midpoint
pts_edge = zeros(length(edges(:,1)),3);
for i = 1:length(edges(:,1))
    pts_edge(i,:) = mean(pts_verts_scaled(edges(i,:),:));
end
% First ghost
pts_ghost1 = pts_edge(1,:);
% Determine position of ghost
dr_ghost1 = bsxfun(@minus,pts_edge,pts_ghost1);
dr_ghost1 = sqrt(sum(dr_ghost1.^2,2));
[~,indx_ghost2] = max(dr_ghost1);
% Second ghost
pts_ghost2 = pts_edge(indx_ghost2,:);
% Place ghost
pts_edge_ghost = [pts_ghost1; pts_ghost2];
% % Define scaling
% for i = 1:length(pts_edge_ghost(:,1))
%     pts_tmp = pts_edge_ghost(i,:);
%     nv_tmp = pts_tmp/norm(pts_tmp);
%     pts_tmp = pts_tmp + nv_tmp*(0.15*(r_insph+sigma_ghost/2));
%     pts_edge_ghost(i,:) = pts_tmp;
% end
% Scale
pts_edge_ghost = ghost_scale_edge*pts_edge_ghost;
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
%%% Build polymer %%%
%%%%%%%%%%%%%%%%%%%%%

% PARAMETER -- Bond location
bond_location = 'corner';

% Define bond location
if strcmp(bond_location,'edge') == 1
    pts_ghost = pts_edge_ghost;
elseif strcmp(bond_location,'face') == 1
    pts_ghost = pts_face_ghost;
elseif strcmp(bond_location,'corner') == 1
    pts_ghost = pts_corner_ghost;
end

% Define rotation - align along z-axis
nv_ghost = pts_ghost(2,:) - pts_ghost(1,:); % dr of line joining 2 ghosts on the same shape centered at origin
nv_ghost = nv_ghost/norm(nv_ghost);

% Define rotation matrix - multiplier for aligning along z-axis
rtn_mtx = rot_az(nv_ghost,1);

% Rotate particle and ghost to final position
pts_verts_final = pts_verts_scaled*rtn_mtx;
pts_ghost_final = pts_ghost*rtn_mtx;

% Define range of growth
if strcmp(shape_name,'Tetrahedron') == 1
    sigma_growth = range(pts_ghost_final(:,3)); % z-axis positions of all the aligned ghost particles
else
    sigma_growth = range(pts_ghost_final(:,3));
end

% Initialize
pts_chain = [];
type_chain = [];
body_chain = [];
q_chain = [];
bond_chain = [];
type_bond_chain = [];

% Type defintion
type_monomer = 0;
type_ghost = 1;

% Buffer factor 
growth_buffer = 1.05;

% Index counter
indx_counter = 0;
% Grow chain
for i = 1:N
    % Base quaternion
    if strcmp(shape_name,'Tetrahedron') == 1
        if mod(i,2) == 0
            if strcmp(bond_location,'face') == 1
                % Rotation
                angle = pi;
                nv_shape_rot = [0 0 1];
                nv_shape_rot = nv_shape_rot/norm(nv_shape_rot);
                q_base = [cos(angle/2) nv_shape_rot*sin(angle/2)];
                sigma_growth = 1.3*range(pts_ghost_final(:,3));
                % Rotate ghost
                pts_ghost_tmp = pts_ghost_final*rot_quat(q_base);
                xshift = pts_ghost_tmp(1,1:2) - pts_ghost_final(2,1:2);
            else
                q_base = [1 0 0 0];
                xshift = 0;
            end            
        else
            nv_growth = [0 0 1];
            nv_growth = nv_growth/norm(nv_growth);
            q_base = [1 0 0 0];
        end
    else
        nv_growth = [0 0 1];
        nv_growth = nv_growth/norm(nv_growth);
        q_base = [1 0 0 0];
    end
    % Define base position
    if i == 1
        pts_base = [0 0 0];
    else
        % Find last monomer
        indx_monomer = find(type_chain == type_monomer);
        indx_monomer = max(indx_monomer);
        pts_base = pts_chain(indx_monomer,:);
        pts_base = pts_base + nv_growth*(sigma_growth + growth_buffer*sigma_ghost);
        if strcmp(shape_name,'Tetrahedron') == 1
            if i > 1
                if mod(i,2) == 0
                    pts_base(1:2) = pts_base(1:2) - xshift;
                else
                    pts_base(1:2) = pts_base(1:2) + xshift;
                end
            end
        end
            
    end    
    % Shift ghost and core to location
    pts_core_add = pts_base;
    pts_ghost_add = bsxfun(@plus,pts_ghost_final*rot_quat(q_base),pts_base);
    % Define index (for body tag)
    indx_body = indx_counter;
    % Update for first ghost    
    indx_counter = indx_counter + 1;
    indx_ghost1 = indx_counter;
    % Update for second ghost
    indx_counter = indx_counter + 1;
    indx_ghost2 = indx_counter;
    % Update for next core
    indx_counter = indx_counter + 1;
    % Create add structures
    pts_add = [pts_core_add; pts_ghost_add];
    type_add = [type_monomer; type_ghost; type_ghost];
    body_add = [indx_body; indx_body; indx_body];
    q_add = [q_base; q_base; q_base];  
    % Bond
    if i > 1
        bond_add = [(indx_body-1) indx_ghost1];
        bond_chain = [bond_chain; bond_add];
        type_bond_chain = [type_bond_chain; 0];
    end
    % Store
    pts_chain = [pts_chain; pts_add];
    type_chain = [type_chain; type_add];
    body_chain = [body_chain; body_add];
    q_chain = [q_chain; q_add];
end

% Center
pts_chain = bsxfun(@minus,pts_chain,mean(pts_chain));

% Box
Lx = max(range(pts_chain(:,1:2)));
Ly = max(range(pts_chain(:,1:2)));
Lz = max(range(pts_chain(:,3)));

% Buffer
Lx = Lx + (2*r_circum);
Ly = Ly + (2*r_circum);
Lz = Lz + (2*r_circum);

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Volume calculation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total volume
V_total = v_shape*N + (length(pts_chain(:,1))-N)*v_ghost/2;

% Target monomer density
rho_mono = 0.85;
rho_mono_sigma = rho_mono/(2*r_avg)^3;

% Average N
Navg = N*(v_shape*N/V_total) + (length(pts_chain(:,1))-N)*((length(pts_chain(:,1))-N)*v_ghost/2/V_total);

% Box size
L = (N/(rho_mono_sigma))^(1/3);

% Volume fraction
phi = V_total/L^3

% return

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
%%% Writing to file %%%
%%%%%%%%%%%%%%%%%%%%%%%

% Write initial configuration
file_name = 'init.txt';
fileID = fopen(file_name,'w');
for i = 1:length(pts_chain(:,1))
    tmp_write = [type_chain(i) body_chain(i) pts_chain(i,:) q_chain(i,:) (2*r_insph) (2*r_avg) sigma_ghost Lx Ly Lz];
    fprintf(fileID,'%i %i %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e\n',tmp_write);
end
fclose(fileID);

% Write initial configuration
file_name = 'init_bond.txt';
fileID = fopen(file_name,'w');
for i = 1:length(bond_chain(:,1))
    tmp_write = [type_bond_chain(i) bond_chain(i,:)];
    fprintf(fileID,'%i %i %i\n',tmp_write);
end
fclose(fileID);

% Write shape vertices
file_name = 'verts.txt';
fileID = fopen(file_name,'w');
for i = 1:length(pts_verts_final(:,1))
    tmp_write = [pts_verts_final(i,:) v_shape I_tensor(1,1) I_tensor(2,2) I_tensor(3,3)];
    fprintf(fileID,'%12e %12e %12e %12e %12e %12e %12e\n',tmp_write);
end
fclose(fileID);

% Write ghost (rigid body)
file_name = 'verts_rigid.txt';
fileID = fopen(file_name,'w');
for i = 1:length(pts_ghost_final(:,1))
    tmp_write = [type_ghost pts_ghost_final(i,:) 1 0 0 0];
    fprintf(fileID,'%i %12e %12e %12e %12e %12e %12e %12e\n',tmp_write); 
end
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

return

%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting check %%%
%%%%%%%%%%%%%%%%%%%%%%

% Create ghost sphere
[xs,ys,zs] = sphere(50);

% Ghost
xs_ghost = (sigma_ghost/2)*xs;
ys_ghost = (sigma_ghost/2)*ys;
zs_ghost = (sigma_ghost/2)*zs;
ghost_sph = [xs_ghost(:) ys_ghost(:) zs_ghost(:)];

text_name = {};
for i = 1:length(pts_verts_scaled(:,1))
    text_name{i} = num2str(i);
end

purple = [107 62 154]/255;
blue = [0 64 128]/255;
red = [180 7 5]/255;
green = [0.4660, 0.6740, 0.1880];
gray = [150 150 150]/255;

figure(1)
face_center = face_center*rtn_mtx;
hold on;
scatter3(pts_verts_final(:,1),pts_verts_final(:,2),pts_verts_final(:,3),100,'filled')
text(pts_verts_final(:,1),pts_verts_final(:,2),pts_verts_final(:,3),text_name,'FontSize',30)
% plot(alphaShape(pts_verts_final,1E10),'EdgeColor','none','FaceAlpha',0.5,'FaceColor',blue)   
scatter3(face_center(:,1),face_center(:,2),face_center(:,3),500,'filled')
for i = 1:length(pts_ghost_final(:,1))
    pts_tmp = bsxfun(@plus,ghost_sph,pts_ghost_final(i,:));
    shp_tmp = alphaShape(pts_tmp,1E10);
    plot(shp_tmp,'EdgeColor','none','FaceColor',purple)
end
% % Only plot ghost
% if strcmp(bond_location,'face') == 1
%     pts_face_ghost_rot = pts_face_ghost*rtn_mtx;
%     for i = 1:length(pts_face_ghost(:,1))
%         pts_tmp = bsxfun(@plus,ghost_sph,pts_face_ghost_rot(i,:));
%         shp_tmp = alphaShape(pts_tmp,1E10);
%         plot(shp_tmp,'EdgeColor','none','FaceColor',purple)
%     end
% elseif strcmp(bond_location,'corner') == 1
%     pts_corner_ghost_rot = pts_corner_ghost*rtn_mtx;
%     for i = 1:length(pts_corner_ghost(:,1))
%         pts_tmp = bsxfun(@plus,ghost_sph,pts_corner_ghost_rot(i,:));
%         shp_tmp = alphaShape(pts_tmp,1E10);
%         plot(shp_tmp,'EdgeColor','none','FaceColor',red)
%     end
% elseif strcmp(bond_location,'edge') == 1
%     pts_edge_ghost_rot = pts_edge_ghost*rtn_mtx;
%     for i = 1:length(pts_edge_ghost(:,1))
%         pts_tmp = bsxfun(@plus,ghost_sph,pts_edge_ghost_rot(i,:));
%         shp_tmp = alphaShape(pts_tmp,1E10);
%         plot(shp_tmp,'EdgeColor','none','FaceColor',green)
%     end
% end
for i = 1:length(edges(:,1))
    ktmp = edges(i,:);   
    plot3(pts_verts_final(ktmp,1),pts_verts_final(ktmp,2),pts_verts_final(ktmp,3),'k','Linewidth',3.5)
end
% for i = 1:length(faces)
%     ktmp = faces{i};
%     ktmp = [ktmp ktmp(1)];
%     plot3(pts_verts_scaled(ktmp,1),pts_verts_scaled(ktmp,2),pts_verts_scaled(ktmp,3),'k','Linewidth',3.5)
% end
% camlight
% lighting flat
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

figure(2)
hold on;
for i = 1:length(pts_chain(:,1))
    pts_tmp = pts_chain(i,:);
    if type_chain(i) == type_monomer
        pts_verts_tmp = pts_verts_final*rot_quat(q_chain(i,:));
        color = blue;
    elseif type_chain(i) == type_ghost
        pts_verts_tmp = ghost_sph;
        color = purple;
    end
    pts_verts_tmp = bsxfun(@plus,pts_verts_tmp,pts_tmp);
    shp_tmp = alphaShape(pts_verts_tmp,1E10);
    plot(shp_tmp,'EdgeColor','none','FaceColor',color,'FaceAlpha',0.2);
end
for i = 1:length(bond_chain(:,1))
    ktmp = bond_chain(i,:)+1;
    plot3(pts_chain(ktmp,1),pts_chain(ktmp,2),pts_chain(ktmp,3),'k','Linewidth',3.5)
end
camlight
lighting flat
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%