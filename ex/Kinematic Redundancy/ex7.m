clear all
clc

syms q1 q2 q3

% Direct kinematics of 3R planar robot
p = [cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
     sin(q1)+sin(q1+q2)+sin(q1+q2+q3);]

% Jacobian 3R
J = jacobian(p,[q1,q2,q3])
J_subs = double(subs(J, [q1,q2,q3],[pi/4,0,pi/4,]))

J_1 = [1 1 1;];
J_extended = [J_1; J_subs;]
rank_J_extended = rank(J_extended)

v = [0; 2; -1;];
J_extended_pinv = pinv(J_extended);
q_dot_ps = J_extended_pinv*v

e_ps = v - J_extended*q_dot_ps

I = [1 0 0;
     0 1 0;
     0 0 1;];
J_1_pinv = pinv(J_1);
v_1 = [0;];
v_2 = [2; -1;];

P_1 = I - J_1_pinv*J_1
q_dot_tp = J_1_pinv*v_1 + pinv(J_subs*P_1)*(v_2 - J_subs*J_1_pinv*v_1)

e_tp = v - J_extended*q_dot_tp

norm_e_ps = norm(e_ps)
norm_e_tp = norm(e_tp)


