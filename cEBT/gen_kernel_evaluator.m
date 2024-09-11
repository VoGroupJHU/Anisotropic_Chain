function [Fk] = gen_kernel_evaluator(pts_coord,shape_grid)

% Generate shape
pts_poly = gen_poly(pts_coord,shape_grid);

% Generate kernel
grid_step = 0.025;
[pos_theo,~,~,theta_space_poly,phi_space_poly] = parameterize_shape(pts_poly,grid_step);
kernel = sqrt(pos_theo{1}.^2 + pos_theo{2}.^2 + pos_theo{3}.^2);
kernel = kernel/min(kernel(:));

% Kernel grid
[tt,pp] = meshgrid(theta_space_poly,phi_space_poly);

% Replicate image array
rep_img = [-1 0 1];

% Create periodic images
kernel_img = [];
% Loop for theta
for i = 1:length(rep_img)
    % Image
    theta_img = rep_img(i);
    % Replication
    tt_rep = tt + theta_img*(2*pi);
    % Loop for phi
    for j = 1:length(rep_img)
        phi_img = rep_img(j);
        pp_rep = pp + phi_img*pi;
        % Storing
        if j == 1
            tmp = [tt_rep(:) -pp_rep(:) kernel(:)];
        elseif j == 2
            tmp = [tt_rep(:) pp_rep(:) kernel(:)];
        elseif j == 3
            tmp = [tt_rep(:) -pp_rep(:) kernel(:)];
        end
        kernel_img = [kernel_img; tmp];
    end
end

% Remove overlap points
kernel_img = roundn(kernel_img,-4);
[~,indx_img] = unique(kernel_img(:,1:2),'stable','rows');
kernel_img = kernel_img(indx_img,:);

% Define step size
dtheta = theta_space_poly(2)-theta_space_poly(1);
dphi = phi_space_poly(2)-phi_space_poly(1);

% Define bracket range to use
%%% NOTE: Take n steps back and forward in each direction %%%
%%% This ensure no NaN upon extrapolation %%%
n_step = 3;
% phi range
phi_low = -pi/2 - n_step*dphi;
phi_high = pi/2 + n_step*dphi;
% theta range
theta_low = -pi - n_step*dtheta;
theta_high = pi + n_step*dtheta;
% Find inde
indx = find( (kernel_img(:,1) >= theta_low) & (kernel_img(:,1) <= theta_high) & (kernel_img(:,2) >= phi_low) & (kernel_img(:,2) <= phi_high));

% Grab relevant point
kernel_table = kernel_img(indx,:);
kernel_table = sortrows(kernel_table, [2 1]);

% Convert to gridded data
Ntheta = length(unique(kernel_table(:,1)));
Nphi = length(unique(kernel_table(:,2)));
tt_table = reshape(kernel_table(:,1),Ntheta,Nphi);
pp_table = reshape(kernel_table(:,2),Ntheta,Nphi);
kk_table = reshape(kernel_table(:,3),Ntheta,Nphi);

% Create kernel evaluator
Fk = griddedInterpolant(tt_table,pp_table,kk_table);