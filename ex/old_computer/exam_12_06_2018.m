clear all
clc

syms q1_dot q2_dot q3_dot q4_dot L1 L2 L3 L4 dc1 dc2 dc3 dc4 q1 q2 q3 q4

%% link 2
v2 = [-L1*sin(q1)*q1_dot-dc2*sin(q2)*q2_dot;
      L1*cos(q1)*q1_dot+dc2*cos(q2)*q2_dot;]
v2x = v2(1,1)
v2y = v2(2)

v2_quad = simplify(((v2x)^2+(v2y)^2), Steps=100)

%% link 3

v3 = [-(L1*sin(q1)*q1_dot+L2*sin(q2)*q2_dot+dc3*sin(q3)*q3_dot);
      L1*cos(q1)*q1_dot+L2*cos(q2)*q2_dot+dc3*cos(q3)*q3_dot;] 

v3x = v3(1,1)
v3y = v3(2)

v3_quad = simplify(((v3x)^2+(v3y)^2), Steps=100)