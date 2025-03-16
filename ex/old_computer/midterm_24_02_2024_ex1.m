clear all
clc

syms q1 q2 q3 l 

p = [l*cos(q1)+l*cos(q1+q2)+l*cos(q1+q2+q3);
     l*sin(q1)+l*sin(q1+q2)+l*sin(q1+q2+q3);]

J = jacobian(p,[q1,q2,q3])

J_subs = double(subs(J, [q1,q2,q3,l],[0,pi/2,-pi/2,1]))

J_pseudoinverse_subs = pinv(J_subs)