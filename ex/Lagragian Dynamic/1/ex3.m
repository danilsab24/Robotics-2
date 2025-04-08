syms q1 q2 q3 l2 real    % joint coordinates
syms dq1 dq2 dq3 real % joint velocities
syms d I2 real

% Joint 2
x2 = d*cos(q2);
y2 = q1 + d*sin(q2);
z2 = 0;  

% Velocity components (vx2, vy2, vz2)
vx2 = diff(x2, q1)*dq1 + diff(x2, q2)*dq2 + diff(x2, q3)*dq3;
vy2 = diff(y2, q1)*dq1 + diff(y2, q2)*dq2 + diff(y2, q3)*dq3;
vz2 = diff(z2, q1)*dq1 + diff(z2, q2)*dq2 + diff(z2, q3)*dq3;

disp('Joint 2:');
disp('p2 = [x2; y2; z2] =');
disp([x2; y2; z2]);
disp('Velocity (vx2, vy2, vz2) =');
disp([vx2; vy2; vz2]);

vect1 = [vx2 vy2 vz2];
vect2 = [vx2; vy2; vz2];
pro   = vect1 * vect2;

% Joint 3
x3 = l2*cos(q2) - q3*sin(q2);
y3 = q1 + l2*sin(q2) + q3*cos(q2);
z3 = 0;  

% Velocity components (vx3, vy3, vz3)
vx3 = diff(x3, q1)*dq1 + diff(x3, q2)*dq2 + diff(x3, q3)*dq3;
vy3 = diff(y3, q1)*dq1 + diff(y3, q2)*dq2 + diff(y3, q3)*dq3;
vz3 = diff(z3, q1)*dq1 + diff(z3, q2)*dq2 + diff(z3, q3)*dq3;

disp('Joint 3:');
disp('p3 = [x3; y3; z3] =');
disp([x3; y3; z3]);
disp('Velocity (vx3, vy3, vz3) =');
disp([vx3; vy3; vz3]);

vect3 = [vx3 vy3 vz3];
vect4 = [vx3; vy3; vz3];
pr3 = simplify(vect3 * vect4, 'Steps',100);

% INERTIA MATRIX M(q)
syms m1 m2 m3 I3 real

% Total kinetic energy:
%  T = (1/2)*m1*(dq1^2)
%      + (1/2)*m2*pro    + (1/2)*I2*(dq2^2)
%      + (1/2)*m3*pr3    + (1/2)*I3*(dq2^2)
%
%  'pro' and 'pr3' are the squared speeds of link2's and link3's CoM
%  computed above.

T = 0.5*m1*(dq1^2) ...
  + 0.5*m2*pro + 0.5*I2*(dq2^2) ...
  + 0.5*m3*pr3 + 0.5*I3*(dq2^2);

% Compute the 3×3 inertia matrix via second derivatives of T w.r.t. dq1,dq2,dq3
M11 = diff(diff(T, dq1), dq1);
M12 = diff(diff(T, dq1), dq2);
M13 = diff(diff(T, dq1), dq3);

M21 = diff(diff(T, dq2), dq1);
M22 = diff(diff(T, dq2), dq2);
M23 = diff(diff(T, dq2), dq3);

M31 = diff(diff(T, dq3), dq1);
M32 = diff(diff(T, dq3), dq2);
M33 = diff(diff(T, dq3), dq3);

M = [M11, M12, M13;
     M21, M22, M23;
     M31, M32, M33];

M = simplify(M);

disp('The inertia matrix M(q) is:');
disp(M);

%% ------------------------------
% CORIOLIS/CENTRIFUGAL VECTOR c(q,dq)
% using the Christoffel symbols for a 3-DOF system
q  = [q1; q2; q3];
dq = [dq1; dq2; dq3];

% Inizializza la matrice C(q, dq) (3x3) a zeri simbolici:
C = sym(zeros(3,3));

% Formula dei coefficienti di Christoffel (Robotics):
%   C(i,j) = 1/2 * SUM_k [ dM(i,j)/dq_k + dM(i,k)/dq_j - dM(k,j)/dq_i ] * dq_k
% Poi il vettore c(q,dq) = C(q,dq) * dq
for i = 1:3
    for j = 1:3
        % Costruisco ogni elemento C(i,j):
        tmp = 0;
        for k = 1:3
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