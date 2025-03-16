clear all
clc

J = [3     1  2;
     1.5  0.5 1;];
J_pinv = pinv(J);

x_dot = [2; 1;]

q_dot = J_pinv*x_dot

q_dot_star_norm = norm(q_dot)

q_dot_prime = [0; 0; 1;];
q_dot_prime_norm = norm(q_dot_prime)