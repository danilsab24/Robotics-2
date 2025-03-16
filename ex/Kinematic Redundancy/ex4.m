clear all
clc

syms q1 q2 q3

% Direct kinematics of 3R planar robot
p = [cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
     sin(q1)+sin(q1+q2)+sin(q1+q2+q3);]

% Jacobian 3R
J = jacobian(p,[q1,q2,q3])
J_subs = double(subs(J, [q1,q2,q3],[pi/2,pi/3,(-2/3)*pi]))

% Objective function
H_q = (sin(q2))^2 + (sin(q3))^2
grad_H_q = [0;
            diff(H_q,q2);
            diff(H_q,q3);]
grad_H_q_subs = double(subs(grad_H_q,[q1,q2,q3],[pi/2,pi/3,(-2/3)*pi]))

% Reduced Gradient methods
J_a = J_subs(:,[1,3])
J_b = J_subs(:, 2)
J_a_inv = inv(J_a)
v = [1; -sqrt(3)]

J_a_inv_algo = [J_a_inv; zeros(1, size(J_a_inv,2))];
pro = -J_a_inv*J_b
pro_J_a_J_b = [pro; ones(1, size(pro,2))]
pro_2 = pro.'
pro_2_J_a_J_b = [pro_2 ones(1, size(pro_2,1))]
final = pro_J_a_J_b*pro_2_J_a_J_b*grad_H_q_subs
q_dot = J_a_inv_algo*v + final

% !! NOTA !! i valori rappresentano (q1 q3 q2) 
% per questo che risultano diversi dalle soluzioni de prof :) 

%% Extended Jacobian and Task Augmentation
J_e = [-(sin(q1)+sin(q1+q2)+sin(q1+q2+q3))      -sin(q1+q2)-sin(q1+q2+q3)   -sin(q1+q2+q3);
        cos(q1)+cos(q1+q2)+cos(q1+q2+q3)        cos(q1+q2)+cos(q1+q2+q3)    cos(q1+q2+q3);
        2*(cos(q1)+cos(q1+q2))                  -2*sin(q2)-3*cos(q1+q2)            0;]

J_e_subs = double(subs(J_e,[q1,q2,q3],[pi/2,pi/3,(-2/3)*pi]))

det_J_e = det(J_e_subs)