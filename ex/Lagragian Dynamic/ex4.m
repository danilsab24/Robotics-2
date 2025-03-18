clear all
clc

syms q1 q2 q3 q4
syms dq1 dq2 dq3 dq4 real

%% Kinetic Energy joint 1
% Position p3 = (x, y, z)
x1 = q1;
y1 = 0;

% Velocity components (vx, vy)
vx1 = diff(x1, q1)*dq1 + diff(x1, q2)*dq2 + diff(x1, q3)*dq3+ diff(x1, q4)*dq4;
vy1 = diff(y1, q1)*dq1 + diff(y1, q2)*dq2 + diff(y1, q3)*dq3+ diff(x1, q4)*dq4;

% Kinetic Energy
T1 = [vx1 vy1]*[vx1; vy1;]

%% Kinetic Energy joint 2
% Position p3 = (x, y, z)
x2 = q1;
y2 = q2;

% Velocity components (vx, vy)
vx2 = diff(x2, q1)*dq1 + diff(x2, q2)*dq2 + diff(x2, q3)*dq3+ diff(x2, q4)*dq4;
vy2 = diff(y2, q1)*dq1 + diff(y2, q2)*dq2 + diff(y2, q3)*dq3+ diff(x2, q4)*dq4;

% Kinetic Energy

%% Kinetic Energy joint 3
% Position p3 = (x, y, z)
x3 = q1+q3;
y3 = q2;

% Velocity components (vx, vy)
vx3 = diff(x3, q1)*dq1 + diff(x3, q2)*dq2 + diff(x3, q3)*dq3+ diff(x4, q4)*dq4;
vy3 = diff(y3, q1)*dq1 + diff(y3, q2)*dq2 + diff(y3, q3)*dq3+ diff(x4, q4)*dq4;

% Kinetic Energy

%% Kinetic Energy joint 4
% Position p3 = (x, y, z)
x4 = q1+q3;
y4 = q2+q4;

% Velocity components (vx, vy)
vx4 = diff(x4, q1)*dq1 + diff(x4, q2)*dq2 + diff(x4, q3)*dq3+ diff(x4, q4)*dq4;
vy4 = diff(y4, q1)*dq1 + diff(y4, q2)*dq2 + diff(y4, q3)*dq3+ diff(x4, q4)*dq4;

% Kinetic Energy


