clear all 
clc

J = [3      1       2;
     1.5    0.5     1;]

x_dot = [2; 1;]

J_pseudoINV = pinv(J)

q_dot = J_pseudoINV*x_dot


