clear all
clc

syms q1 q2 q3 q4 l

p = [l*cos(q1)+l*cos(q1+q2)+l*cos(q1+q2+q3)+l*cos(q1+q2+q3+q4);
     l*sin(q1)+l*sin(q1+q2)+l*sin(q1+q2+q3)+l*sin(q1+q2+q3+q4);
     q1+q2+q3+q4];

J = jacobian(p,[q1,q2,q3,q4]);
disp('Expression of Jacobian: ');
disp(J)

J_subs = double(subs(J,[q1,q2,q3,q4,l],[0,0,pi/2,0,0.5]));
disp('Jacobian = ');
disp(J_subs)

J_subs_pinv = pinv(J_subs);
disp('Pseudoinverse Weighted Jacobian = ');
disp(J_subs_pinv)

% Number of joints
n = 4;

% Joint range: q_i in [-2, 2] for each i
q_min = -2*ones(n,1);
q_max =  2*ones(n,1);

% Mid-range q_bar (which is 0 for all joints in this case)
q_bar = (q_min + q_max)/2;  % = [0;0;0;0]

% H_range(q) = 1/(2*n) * SUM[ ((q_i - q_bar_i)/(q_max_i-q_min_i))^2 ]
% Here, for each joint, (q_bar_i=0) and (q_max_i-q_min_i=4),
% so the term becomes (q_i/4)^2

H_range = (1/(2*n)) * ( (q1/4)^2 ...
                      + (q2/4)^2 ...
                      + (q3/4)^2 ...
                      + (q4/4)^2 );

% Gradient of H_range w.r.t. [q1, q2, q3, q4]
grad_Hrange = [ diff(H_range, q1);
                diff(H_range, q2);
                diff(H_range, q3);
                diff(H_range, q4) ];

% Example joint configuration at which we want to evaluate
q_val = [ 0, 0, pi/2, 0];  % you can pick any feasible point in [-2,2]

% Evaluate H_range and its gradient numerically at q_val
H_val     = double(subs(H_range,     [q1,q2,q3,q4], q_val));
gradH_val = double(subs(grad_Hrange, [q1,q2,q3,q4], q_val));

% Display results
disp('H_range(q_val) =');
disp(H_val);

disp('Gradient of H_range(q_val) =');
disp(gradH_val);

% Compute q dot
v = [1; 0; 0.5];
I = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];
q_r = J_subs_pinv*v
q_n = - (I - J_subs_pinv*J_subs)*gradH_val;
q_dot =  q_r + q_n;
disp('Value of q r = ');
disp(q_r)
disp('Value of q n = ');
disp(q_n)
disp('Value of q dot = ');
disp(q_dot)

% In case of joint velocity out of bound
q_max = 0.5;
k = (q_max-q_n(1))/q_r(1);
disp('Value of k =');
disp(k);

q_dot_scale = q_r*k + q_n;
disp('Value of q dot scaled = ');
disp(q_dot_scale)