clear all
clc

%% 1) Declare symbolic variables

syms q1 q2 q3 real
syms dq1 dq2 dq3 real
syms m1 m2 m3 I3 dC3 real

%% 2) Define the (x,y) positions of each mass, consistent with the figure:

x1 = 0;               y1 = q1;
x2 = q2;               y2 = q1;
x3 = q2 + dC3*cos(q3); y3 = q1 + dC3*sin(q3);

%% 3) Compute the velocities of each mass (vx, vy) by taking time-derivatives

% Velocity of m1
v1x = diff(x1,q1)*dq1 + diff(x1,q2)*dq2 + diff(x1,q3)*dq3;
v1y = diff(y1,q1)*dq1 + diff(y1,q2)*dq2 + diff(y1,q3)*dq3;

% Velocity of m2
v2x = diff(x2,q1)*dq1 + diff(x2,q2)*dq2 + diff(x2,q3)*dq3;
v2y = diff(y2,q1)*dq1 + diff(y2,q2)*dq2 + diff(y2,q3)*dq3;

% Velocity of m3
v3x = diff(x3,q1)*dq1 + diff(x3,q2)*dq2 + diff(x3,q3)*dq3;
v3y = diff(y3,q1)*dq1 + diff(y3,q2)*dq2 + diff(y3,q3)*dq3;

%% 4) Translational kinetic energy for each mass: (1/2)*m*(vx^2 + vy^2)
T1_trans = 0.5 * m1 * (v1x^2 + v1y^2);
T2_trans = 0.5 * m2 * (v2x^2 + v2y^2);
T3_trans = 0.5 * m3 * (v3x^2 + v3y^2);

%% 5) Rotational kinetic energy for link 3

T3_rot = 0.5 * I3 * dq3^2;

%% 6) Total kinetic energy T
T = T1_trans + T2_trans + (T3_trans + T3_rot);

%% 7) Build the inertia matrix M(q) by second derivatives of T wrt dq1, dq2, dq3
%    M(i,j) = d^2 T / (d dq_i d dq_j).
M11 = diff(diff(T, dq1), dq1);
M12 = diff(diff(T, dq1), dq2);
M13 = diff(diff(T, dq1), dq3);

M21 = diff(diff(T, dq2), dq1);
M22 = diff(diff(T, dq2), dq2);
M23 = diff(diff(T, dq2), dq3);

M31 = diff(diff(T, dq3), dq1);
M32 = diff(diff(T, dq3), dq2);
M33 = diff(diff(T, dq3), dq3);

% Combine into a matrix
M = [M11, M12, M13;
     M21, M22, M23;
     M31, M32, M33];

% Simplify the symbolic expressions
M = simplify(M);

%% 8) Display the final inertia matrix
disp('The inertia matrix M(q) is:');
disp(M);
