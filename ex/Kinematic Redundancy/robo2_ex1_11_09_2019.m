clear all
clc

syms q1 q2 q3 L

p = [L*cos(q1)+L*cos(q1+q2)+L*cos(q1+q2+q3);
     L*sin(q1)+L*sin(q1+q2)+L*sin(q1+q2+q3);]

J = simplify(jacobian(p,[q1,q2,q3]),Steps=100)