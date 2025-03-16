clear all
clc

% exercise EXTENDED JACOBIAN PSEUDO INVERSE METHOD DLS TASK PROIORITY

syms q1 q2 q3

p = [cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
     sin(q1)+sin(q1+q2)+sin(q1+q2+q3);]

J = jacobian(p,[q1,q2,q3])

J_subs = double(subs(J,[q1,q2,q3],[pi/4,0,0]))

J_pinv_subs = pinv(J_subs)