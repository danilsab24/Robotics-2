clear all
clc

syms q1 q2 q3 q4

p_e = [cos(q1)+cos(q2)+cos(q3)+cos(q4);
       sin(q1)+sin(q2)+sin(q3)+sin(q4);];

J_e = jacobian(p_e,[q1,q2,q3,q4])

J_e_subs = double(subs(J_e,[q1,q2,q3,q4],[0,pi/6,-pi/3,-pi/3]))

p_t = [cos(q1)+cos(q2);
       sin(q1)+sin(q2);]

J_t = jacobian(p_t,[q1,q2,q3,q4]);

J_t_subs = double(subs(J_t,[q1,q2,q3,q4],[0,pi/6,-pi/3,-pi/3]))

rank_J_e = rank(J_e_subs)

rank_J_t = rank(J_t_subs)

J_complete = [J_e_subs; J_t_subs;]
rank_J_complete = rank(J_complete)

J_e_pinv = pinv(J_e_subs);
J_t_pinv = pinv(J_t_subs);
J_complete_pinv = pinv(J_complete);

v_e = [0.4330; -0.75;];
v_t = [-0.5; 0.8660;];
v_complete = [v_e;v_t;];

q_a_dot = J_e_pinv*v_e
q_b_dot = J_t_pinv*v_t

e1_e = v_e - (J_e_subs*q_a_dot);
e1_t = v_t - (J_t_subs*q_a_dot)
e1_t_norm = norm(e1_t)

e2_e = v_e - (J_e_subs*q_b_dot)
e2_t = v_t - (J_t_subs*q_b_dot);
e2_t_norm = norm(e2_e)

q_c_dot = J_complete_pinv*v_complete
e3_e = v_complete - (J_complete*q_c_dot)

I = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1;];
P_e = I - J_e_pinv*J_e_subs;
A = J_e_pinv*v_e;
B =  pinv(J_t_subs*P_e);
C = (v_t - J_t_subs*J_e_pinv*v_e);
q_d_dot = A + B*C


P_t = I - J_t_pinv*J_t_subs;
A2 = J_t_pinv*v_t;
B2 =  pinv(J_e_subs*P_t);
C2 = (v_e - J_e_subs*J_t_pinv*v_t);
q_e_dot = A2 + B2*C2


