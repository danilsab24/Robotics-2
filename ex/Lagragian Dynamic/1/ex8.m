clear all
clc

syms q1 q2 q3 l1 l2 l3 m1 m2 m3 d1 d2 d3  I1 I2 I3 g0 real
syms dq1 dq2 dq3 ddq1 ddq2 real

%% Kinetic energy of joint 1
y1 = 0;
x1 = q1;
z1 = 0;

vx1 = diff(x1, q1)*dq1 + diff(x1, q2)*dq2 + diff(x1, q3)*dq3;
vy1 = diff(y1, q1)*dq1 + diff(y1, q2)*dq2 + diff(y1, q3)*dq3;
vz1 = diff(z1, q1)*dq1 + diff(z1, q2)*dq2 + diff(z1, q3)*dq3;

T1 = 0.5*m1*[vx1 vy1 0]*[vx1; vy1; 0;]+0.5*I1*[0 0 vz1]*[0; 0; vz1;];

disp('velocity of joint 1: ');
disp([vx1; vy1; vz1;]);

disp('Kinetic Energy of joint 1: ');
disp(T1);

%% Kinetic energy of joint 2
y2 = d2*sin(q2);
x2 = q1+d2*cos(q2);
z2 = q2;

vx2 = diff(x2, q1)*dq1 + diff(x2, q2)*dq2 + diff(x2, q3)*dq3;
vy2 = diff(y2, q1)*dq1 + diff(y2, q2)*dq2 + diff(y2, q3)*dq3;
vz2 = diff(z2, q1)*dq1 + diff(z2, q2)*dq2 + diff(z2, q3)*dq3;

T2 = 0.5*m2*[vx2 vy2 0;]*[vx2; vy2; 0;] + 0.5*I2*[0 0 vz2;]*[0; 0; vz2;];

disp('velocity of joint 2: ');
disp([vx2; vy2; vz2;]);

disp('Kinetic Energy of joint 2: ');
disp(T2);

T = simplify(T1 + T2, Steps=1000);
disp("Total Kinetic Energy: ");
disp(T);

%% Inertia Matrix
M11 = diff(diff(T, dq1), dq1);
M12 = diff(diff(T, dq1), dq2);

M21 = diff(diff(T, dq2), dq1);
M22 = diff(diff(T, dq2), dq2);

M = [M11, M12;
     M21, M22;];

M = simplify(M);

disp('The inertia matrix M(q) is:');
disp(M);

%% CORIOLIS/CENTRIFUGAL VECTOR c(q,dq)
% using the Christoffel symbols for a 3-DOF system
q  = [q1; q2;];
dq = [dq1; dq2;];

% Inizializza la matrice C(q, dq) (2x2) a zeri simbolici:
C = sym(zeros(2,2));

% Formula dei coefficienti di Christoffel (Robotics):
%   C(i,j) = 1/2 * SUM_k [ dM(i,j)/dq_k + dM(i,k)/dq_j - dM(k,j)/dq_i ] * dq_k
% Poi il vettore c(q,dq) = C(q,dq) * dq
for i = 1:2
    for j = 1:2
        % Costruisco ogni elemento C(i,j):
        tmp = 0;
        for k = 1:2
            tmp = tmp + 0.5 * ( ...
                diff(M(i,j), q(k)) + ...
                diff(M(i,k), q(j)) - ...
                diff(M(k,j), q(i)) ) * dq(k);
        end
        C(i,j) = simplify(tmp);  % semplifico
    end
end

% Finalmente calcolo il vettore c(q, dq):
c_vec = simplify(C * dq);

disp('La matrice di Coriolis C(q,dq) è:');
disp(C);

disp('Il vettore c(q,dq) = C(q,dq)*dq è:');
disp(c_vec);

%% Gravity terms
U1 = m1*g0*x1;
U2 = m2*g0*x2;

U = [U1+U2];
g_q = simplify([diff(U,q1); diff(U,q2);], Steps=100);

disp('Potential Energy: ');
disp(U);
disp('Gravity terms g(q): ');
disp(g_q);

%% Linearization
tau = M*[ddq1; ddq2;] + g_q;

disp('tau: ');
disp(tau);



