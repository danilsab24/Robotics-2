clear all
clc

syms q1 q2 q3 q4 m1 m2 m3 m4
syms dq1 dq2 dq3 dq4 real

%% Kinetic Energy joint 1
% Position p3 = (x, y, z)
x1 = q1;
y1 = 0;

% Velocity components (vx, vy)
vx1 = diff(x1, q1)*dq1 + diff(x1, q2)*dq2 + diff(x1, q3)*dq3+ diff(x1, q4)*dq4;
vy1 = diff(y1, q1)*dq1 + diff(y1, q2)*dq2 + diff(y1, q3)*dq3+ diff(y1, q4)*dq4;

% Kinetic Energy
T1 = 0.5*m1*[vx1 vy1]*[vx1; vy1;];

%% Kinetic Energy joint 2
% Position p3 = (x, y, z)
x2 = q1;
y2 = q2;

% Velocity components (vx, vy)
vx2 = diff(x2, q1)*dq1 + diff(x2, q2)*dq2 + diff(x2, q3)*dq3+ diff(x2, q4)*dq4;
vy2 = diff(y2, q1)*dq1 + diff(y2, q2)*dq2 + diff(y2, q3)*dq3+ diff(y2, q4)*dq4;

% Kinetic Energy
T2 = 0.5*m2*[vx2 vy2]*[vx2; vy2;];

%% Kinetic Energy joint 3
% Position p3 = (x, y, z)
x3 = q1+q3;
y3 = q2;

% Velocity components (vx, vy)
vx3 = diff(x3, q1)*dq1 + diff(x3, q2)*dq2 + diff(x3, q3)*dq3+ diff(x3, q4)*dq4;
vy3 = diff(y3, q1)*dq1 + diff(y3, q2)*dq2 + diff(y3, q3)*dq3+ diff(y3, q4)*dq4;

% Kinetic Energy
T3 = 0.5*m3*[vx3 vy3]*[vx3; vy3;];

%% Kinetic Energy joint 4
% Position p3 = (x, y, z)
x4 = q1+q3;
y4 = q2+q4;

% Velocity components (vx, vy)
vx4 = diff(x4, q1)*dq1 + diff(x4, q2)*dq2 + diff(x4, q3)*dq3+ diff(x4, q4)*dq4;
vy4 = diff(y4, q1)*dq1 + diff(y4, q2)*dq2 + diff(y4, q3)*dq3+ diff(y4, q4)*dq4;

% Kinetic Energy
T4 = 0.5*m4*[vx4 vy4]*[vx4; vy4;];

% Sum of all kinetic energy
T = simplify(T1 + T2 + T3 + T4, Steps=1000)

%% Inertia Matrix

M11 = diff(diff(T, dq1), dq1);
M12 = diff(diff(T, dq1), dq2);
M13 = diff(diff(T, dq1), dq3);
M14 = diff(diff(T, dq1), dq4);

M21 = diff(diff(T, dq2), dq1);
M22 = diff(diff(T, dq2), dq2);
M23 = diff(diff(T, dq2), dq3);
M24 = diff(diff(T, dq2), dq4);

M31 = diff(diff(T, dq3), dq1);
M32 = diff(diff(T, dq3), dq2);
M33 = diff(diff(T, dq3), dq3);
M34 = diff(diff(T, dq3), dq4);

M41 = diff(diff(T, dq4), dq1);
M42 = diff(diff(T, dq4), dq2);
M43 = diff(diff(T, dq4), dq3);
M44 = diff(diff(T, dq4), dq4);

M = [M11, M12, M13, M14;
     M21, M22, M23, M24;
     M31, M32, M33, M34;
     M41, M42, M43, M44;];

M = simplify(M);

disp('The inertia matrix M(q) is:');
disp(M);

%% q_dot minimized Kinetic energy Inertia Weightd PseudoInverse
syms q1 q2 q3 q4 v_xd v_yd
pM = [q1+q3;
     q2+q4;];
JM = jacobian(pM,[q1,q2,q3,q4]);
M_inv = inv(M);
v_d_M = [v_xd; v_yd;];
disp('q_dot minimized Kinetic energy Inertia Weightd PseudoInverse: ');
q_dot_M = simplify(((M_inv*JM.')*inv(JM*M_inv*JM.'))*v_d_M, Steps=100)

%% q_dot minimized norm of the joint velocity PseudoInverse

syms q1 q2 q3 q4 v_xd v_yd
p = [q1+q3;
     q2+q4;];
J = jacobian(p,[q1,q2,q3,q4]);
J_pinv = pinv(J);
v_d = [v_xd; v_yd;];
disp('q_dot minimized norm of the joint velocity PseudoInverse');
q_dot = J_pinv*v_d
