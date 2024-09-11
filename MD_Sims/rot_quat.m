function rtn_mtx = rot_quat(q_tmp)
% Paper: "Elucidating the Interplay Between Entropy-Driven and Patch-Mediated 
% Bonding in Directing Nanoscale Assemblies"
% Authors: Kireeti Akkunuri, Xiangyu Zhang, Thi Vo
% Dated: September 2024


% Input
% q_tmp: quaternion

% Output
% rnt_mtx: rotation matrix

% Define quaternion
q0 = q_tmp(1);
q1 = q_tmp(2);
q2 = q_tmp(3);
q3 = q_tmp(4);
% Define rotation matrix
rtn_mtx = zeros(3,3);
% First row
rtn_mtx(1,1) = 1 - 2*q2^2 - 2*q3^2;
rtn_mtx(1,2) = 2*(q1*q2 + q0*q3);
rtn_mtx(1,3) = 2*(q1*q3 - q0*q2);
% Second row
rtn_mtx(2,1) = 2*(q1*q2 - q0*q3);
rtn_mtx(2,2) = 1 - 2*q1^2 - 2*q3^2;
rtn_mtx(2,3) = 2*(q2*q3 + q0*q1);
% Third row
rtn_mtx(3,1) = 2*(q1*q3 + q0*q2);
rtn_mtx(3,2) = 2*(q2*q3 - q0*q1);
rtn_mtx(3,3) = 1 - 2*q1^2 - 2*q2^2;