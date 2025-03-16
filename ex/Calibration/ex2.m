clear all
clc

syms q1 q2 l1 l2

p_a = [2; 0;]; 
p_b = [0; 2;];
p_c = [1.6925; 0.7425;];
p_d = [1.7218; 0.6718;];

q_a = [0; 0;];
q_b = [pi/2; 0;];
q_c = [pi/4; -pi/4;];
q_d = [0; pi/4;];

p = [l1*cos(q1)+l2*cos(q1+q2);
     l1*sin(q1)+l2*sin(q1+q2);];

p_a_hat = double(subs(p, [q1,q2,l1,l2],[0,0,1,1]))
p_b_hat = double(subs(p, [q1,q2,l1,l2],[pi/2,0,1,1]))
p_c_hat = double(subs(p, [q1,q2,l1,l2],[pi/4,-pi/4,1,1]))
p_d_hat = double(subs(p, [q1,q2,l1,l2],[0,pi/4,1,1]))

delta_p_a = p_a - p_a_hat;
delta_p_b = p_b - p_b_hat;
delta_p_c = p_c - p_c_hat;
delta_p_d = p_d - p_d_hat;

delta_p = [delta_p_a; delta_p_b; delta_p_c; delta_p_d;]

phi = [cos(q1) cos(q1+q2);
       sin(q1) sin(q1+q2);];

phi_a = subs(phi, [q1,q2],[0,0]);
phi_b = subs(phi, [q1,q2],[pi/2,0]);
phi_c = subs(phi, [q1,q2],[pi/4,-pi/4]);
phi_d = subs(phi, [q1,q2],[0,pi/4]);

phi_8x2 = double([phi_a; phi_b; phi_c; phi_d])

phi_8x2_psinv = pinv(phi_8x2);

delta_l = phi_8x2_psinv*delta_p
