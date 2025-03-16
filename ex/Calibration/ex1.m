clear all
clc

syms q1 q2 L1 L2

p = [L1*cos(q1)+L2*cos(q1+q2);
     L1*sin(q1)+L2*sin(q1+q2);]

J_theta = jacobian(p, [q1,q2])
J_a = jacobian(p,[L1,L2])