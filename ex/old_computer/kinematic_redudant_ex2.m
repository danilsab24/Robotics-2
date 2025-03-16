clear all
clc

% Define symbolic variables for the three joints
syms q1 q2 q3

% Define the end-effector position for a planar 3R robot
p = [ cos(q1) + cos(q1+q2) + cos(q1+q2+q3);
      sin(q1) + sin(q1+q2) + sin(q1+q2+q3) ];

% Compute the Jacobian matrix
J = simplify(jacobian(p, [q1, q2, q3]), Steps=10)

% Specify the joint values where we want to evaluate the Jacobian
q_val = [(2/5)*pi, (pi/2), (-pi/4)];

% Evaluate the Jacobian numerically at q_val
J_subs = double(subs(J, [q1, q2, q3], q_val));

% Compute the pseudo-inverse of J_subs
J_psinv_subs = pinv(J_subs);

% Display the numerical Jacobian and its pseudo-inverse
disp('J_subs =');
disp(J_subs);

disp('Pseudo-inverse of J_subs =');
disp(J_psinv_subs);

%% ---------------------------------------------------
% Compute the H_range(q) criterion for the given joints
% ---------------------------------------------------
% Joint bounds:
%   q1 in [-pi/2,  pi/2]
%   q2 in [   0 , 2*pi/3]
%   q3 in [-pi/4,  pi/4]
q_min = [-pi/2;     0;       -pi/4 ];
q_max = [ pi/2;  2*pi/3;     pi/4  ];

% Compute the mid-range q_bar for each joint
q_bar = (q_min + q_max) / 2;

% Number of joints
n = 3;

% Define H_range(q) following the figure's formula:
%  H_range(q) = 1/(2n) * sum( ( (q_i - q_bar_i)/(qMax_i - qMin_i) )^2 )
%  with the given ranges and midrange:
%   - q1 in [-pi/2, pi/2] => midrange = 0 => (q1/pi)^2
%   - q2 in [0, 2pi/3]    => midrange = pi/3 => ((q2 - pi/3)/(2pi/3))^2
%   - q3 in [-pi/4, pi/4] => midrange = 0 => (q3/(pi/2))^2
%
% and overall scaling factor 1/(2*3) = 1/6

H_range = (1/6) * ( (q1/pi)^2 ...
                  + ((q2 - pi/3)/(2*pi/3))^2 ...
                  + (q3/(pi/2))^2 );

% Compute the gradient of H_range w.r.t. [q1, q2, q3]
grad_Hrange = [diff(H_range, q1);
               diff(H_range, q2);
               diff(H_range, q3)];

% Joint configuration at which we want to evaluate:
q_val = [ (2/5)*pi,  pi/2,  -pi/4 ];

% Evaluate H_range and its gradient numerically
H_val      = double(subs(H_range,     [q1,q2,q3], q_val));
gradH_val  = double(subs(grad_Hrange, [q1,q2,q3], q_val));

% Display the results
disp('Value of H_range(q) at the given configuration:');
disp(H_val);

disp('Gradient of H_range(q) at the given configuration (as a 3x1 vector):');
disp(gradH_val);

% Compute q dot
v = [-3; 0];
I = [1 0 0;
     0 1 0;
     0 0 1];
q_r = J_psinv_subs*v
q_n = - (I - J_psinv_subs*J_subs)*gradH_val;
q_dot =  q_r + q_n;
disp('Value of q r = ');
disp(q_r)
disp('Value of q n = ');
disp(q_n)
disp('Value of q dot = ');
disp(q_dot)

q_max = 2;
k = (q_max-q_n(1))/q_r(1);
disp('Value of k =');
disp(k);

q_dot_scale = q_r*k + q_n;
disp('Value of q dot scaled = ');
disp(q_dot_scale)
