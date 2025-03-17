clear all
clc

%% 1) Define symbolic variables
%  q1, q2, q3      -> joint angles
%  l1, l2, l3      -> link lengths (DH parameters)
%  m1, m2, m3      -> link masses
%  rcx_i, rcy_i    -> CoM coordinates in the i-th local frame
%  g0              -> magnitude of gravity acting along x0 (downward)
syms q1 q2 q3 real
syms l1 l2 l3 real
syms m1 m2 m3 real
syms rcx1 rcy1 rcx2 rcy2 rcx3 rcy3 real
syms g0 real

%% 2) Define the gravity vector in the chosen base frame
%  According to the problem statement, the gravity vector is (g0, 0, 0)^T 
%  pointing along x0 downward.
g_0 = [g0; 0; 0];

%% 3) Build the individual DH transformation matrices
%  For a 3R planar manipulator (alpha_i = 0, di = 0):
%    i^-1 A_i(q_i) = [ cos(q_i), -sin(q_i),  0,   l_i*cos(q_i)
%                     sin(q_i),  cos(q_i),  0,   l_i*sin(q_i)
%                     0,         0,         1,   0
%                     0,         0,         0,   1 ]
%
%  We will call these T1, T2, T3.

T1 = [ cos(q1), -sin(q1),  0,  l1*cos(q1);
       sin(q1),  cos(q1),  0,  l1*sin(q1);
       0,        0,        1,  0;
       0,        0,        0,  1 ];

T2 = [ cos(q2), -sin(q2),  0,  l2*cos(q2);
       sin(q2),  cos(q2),  0,  l2*sin(q2);
       0,        0,        1,  0;
       0,        0,        0,  1 ];

T3 = [ cos(q3), -sin(q3),  0,  l3*cos(q3);
       sin(q3),  cos(q3),  0,  l3*sin(q3);
       0,        0,        1,  0;
       0,        0,        0,  1 ];

%% 4) Compute the cumulative transformations from the base frame

T_01 = T1;
T_02 = T1 * T2;
T_03 = T1 * T2 * T3;

%% 5) Define the CoM position vectors in each local link frame
r_c1_local = [rcx1; rcy1; 0; 1];
r_c2_local = [rcx2; rcy2; 0; 1];
r_c3_local = [rcx3; rcy3; 0; 1];

%% 6) Transform CoMs into the base frame

r_c1_base = T_01 * r_c1_local;
r_c2_base = T_02 * r_c2_local;
r_c3_base = T_03 * r_c3_local;

%% 7) Compute potential energy for each link

% Extract only the first three coordinates (x0, y0, z0) from r_c_base:
rc1_xyz = r_c1_base(1:3);
rc2_xyz = r_c2_base(1:3);
rc3_xyz = r_c3_base(1:3);

% Potential energies (symbolic):
U1 = - m1 * (g_0.' * rc1_xyz); 
U2 = - m2 * (g_0.' * rc2_xyz); 
U3 = - m3 * (g_0.' * rc3_xyz);

% Simplify them (optional)
U1_simpl = simplify(U1);
U2_simpl = simplify(U2);
U3_simpl = simplify(U3);

disp('--- Symbolic Potential Energy Expressions ---');
disp('U1 = '); disp(U1_simpl);
disp('U2 = '); disp(U2_simpl);
disp('U3 = '); disp(U3_simpl);

%% 8) Compute g(q)
U = [U1+U2+U3];
g_q = simplify([diff(U,q1); diff(U,q2); diff(U,q3);], Steps=100)
