clear all
clc

syms q1 q2 q3 l

% Direct kinematics of 3R planar robot
p = [l*cos(q1)+l*cos(q1+q2)+l*cos(q1+q2+q3);
     l*sin(q1)+l*sin(q1+q2)+l*sin(q1+q2+q3);]

% Jacobian 3R
J = jacobian(p,[q1,q2,q3])
J_subs = double(subs(J, [q1,q2,q3,l],[0,0,pi/2,0.5]))

J_dot_subs = [-0.8 -0.4 0;
                0    0  0;];
n_q_q_dot = [-0.64; 0;];

J_subs_pinv = pinv(J_subs);
p_dot_dot = [2; 1;];
q_dot = [0.8; 0; -0.8;];

q_ddot = J_subs_pinv*(p_dot_dot -J_dot_subs*q_dot)