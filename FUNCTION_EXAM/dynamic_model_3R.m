clear all
clc

syms q1 q2 q3 q_dot_1 q_dot_2 q_dot_3 I_1 I_2 I_3 real
syms l_1 l_2 l_3 d_1 d_2 d_3 m_1 m_2 m_3 g0 positive

% Joint 1
x1 = d_1*cos(q1);
y1 = d_1*sin(q1);

vx1 = diff(x1,q1)*q_dot_1 + diff(x1,q2)*q_dot_2 + diff(x1,q3)*q_dot_3
vy1 = diff(y1,q1)*q_dot_1 + diff(y1,q2)*q_dot_2 + diff(y1,q3)*q_dot_3

T1 = 0.5*m_1*[vx1 vy1]*[vx1; vy1;] + 0.5*I_1*(q_dot_1^2)

% Joint 2
x2 = l_1*cos(q1) + d_2*cos(q1+q2);
y2 = l_1*sin(q1) + d_2*sin(q1+q2);

vx2 = diff(x2,q1)*q_dot_1 + diff(x2,q2)*q_dot_2 + diff(x2,q3)*q_dot_3
vy2 = diff(y2,q1)*q_dot_1 + diff(y2,q2)*q_dot_2 + diff(y2,q3)*q_dot_3

T2 = 0.5*m_2*[vx2 vy2]*[vx2; vy2;] + 0.5*((q_dot_1+q_dot_2)^2)

% Joint 3
x3 = l_1*cos(q1) + l_2*cos(q1+q2) + d_3*cos(q1+q2+q3);
y3 = l_1*sin(q1) + l_2*sin(q1+q2) + d_3*sin(q1+q2+q3);

vx3 = diff(x3,q1)*q_dot_1 + diff(x3,q2)*q_dot_2 + diff(x3,q3)*q_dot_3
vy3 = diff(y3,q1)*q_dot_1 + diff(y3,q2)*q_dot_2 + diff(y3,q3)*q_dot_3

T3 = 0.5*m_3*[vx3 vy3]*[vx3; vy3;] + 0.5*I_3*((q_dot_1+q_dot_2+q_dot_3)^2)

T = T1+T2+T3

% Inertia Matrix
M11 = diff(diff(T, q_dot_1), q_dot_1);
M12 = diff(diff(T, q_dot_1), q_dot_2);
M13 = diff(diff(T, q_dot_1), q_dot_3);

M21 = diff(diff(T, q_dot_2), q_dot_1);
M22 = diff(diff(T, q_dot_2), q_dot_2);
M23 = diff(diff(T, q_dot_2), q_dot_3);

M31 = diff(diff(T, q_dot_3), q_dot_1);
M32 = diff(diff(T, q_dot_3), q_dot_2);
M33 = diff(diff(T, q_dot_3), q_dot_3);

M = simplify([M11, M12, M13; M21, M22, M23; M31, M32, M33])

% Coriolis and Centrifugal terms
q  = [q1; q2; q3;];
dq = [q_dot_1; q_dot_2; q_dot_3;];
C = sym(zeros(3,3));

for i = 1:3
    for j = 1:3
        tmp = 0;
        for k = 1:3
            tmp = tmp + 0.5 * ( ...
                diff(M(i,j), q(k)) + ...
                diff(M(i,k), q(j)) - ...
                diff(M(k,j), q(i)) ) * dq(k);
        end
        C(i,j) = simplify(tmp);  
    end
end

c_vec = simplify(C * dq);

disp('La matrice di Coriolis C(q,dq) è:');
disp(C);
disp('Il vettore c(q,dq) = C(q,dq)*dq è:');
disp(c_vec);

% Gravity terms
% !! NOTE CHANGE y1 y2 and y3 MUST BE CHANGED !!
U1 = m_1*g0*(y1);
U2 = m-2*g0*(y2);
U3 = m_3*g0*(y3);

U = [U1+U2+U3];
g_q = simplify([diff(U,q1); diff(U,q2); diff(U,q3)])