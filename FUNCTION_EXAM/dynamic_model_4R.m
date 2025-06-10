clear all
clc

syms q1 q2 q3 q4 q_dot_1 q_dot_2 q_dot_3 q_dot_4 I_1 I_2 I_3 I_4 real
syms l_1 l_2 l_3 l_4 d_1 d_2 d_3 d_4 m_1 m_2 m_3 m_4 g0 positive

% Joint 1
x1 = d_1*cos(q1);
y1 = d_1*sin(q1);

vx1 = diff(x1,q1)*q_dot_1 + diff(x1,q2)*q_dot_2 + diff(x1,q3)*q_dot_3 + diff(x1,q4)*q_dot_4; 
vy1 = diff(y1,q1)*q_dot_1 + diff(y1,q2)*q_dot_2 + diff(y1,q3)*q_dot_3 + diff(y1,q4)*q_dot_4;

T1 = 0.5*m_1*[vx1 vy1]*[vx1; vy1;] + 0.5*I_1*(q_dot_1^2)

% Joint 2
x2 = l_1*cos(q1) + d_2*cos(q1+q2);
y2 = l_1*sin(q1) + d_2*sin(q1+q2);

vx2 = diff(x2,q1)*q_dot_1 + diff(x2,q2)*q_dot_2 + diff(x2,q3)*q_dot_3 + diff(x2,q4)*q_dot_4;
vy2 = diff(y2,q1)*q_dot_1 + diff(y2,q2)*q_dot_2 + diff(y2,q3)*q_dot_3 + diff(y2,q4)*q_dot_4;

T2 = 0.5*m_2*[vx2 vy2]*[vx2; vy2;] + 0.5*((q_dot_1+q_dot_2)^2)

% Joint 3
x3 = l_1*cos(q1) + l_2*cos(q1+q2) + d_3*cos(q1+q2+q3);
y3 = l_1*sin(q1) + l_2*sin(q1+q2) + d_3*sin(q1+q2+q3); 

vx3 = diff(x3,q1)*q_dot_1 + diff(x3,q2)*q_dot_2 + diff(x3,q3)*q_dot_3 + diff(x3,q4)*q_dot_4;
vy3 = diff(y3,q1)*q_dot_1 + diff(y3,q2)*q_dot_2 + diff(y3,q3)*q_dot_3 + diff(y3,q4)*q_dot_4;

T3 = 0.5*m_3*[vx3 vy3]*[vx3; vy3;] + 0.5*I_3*((q_dot_1+q_dot_2+q_dot_3)^2)

% Joint 4
x4 = l_1*cos(q1) + l_2*cos(q1+q2) + l_3*cos(q1+q2+q3) + d_4*cos(q1+q2+q3+q4);
y4 = l_1*sin(q1) + l_2*sin(q1+q2) + l_3*sin(q1+q2+q3) + d_4*sin(q1+q2+q3+q4); 

vx4 = diff(x4,q1)*q_dot_1 + diff(x4,q2)*q_dot_2 + diff(x4,q3)*q_dot_3 + diff(x4,q4)*q_dot_4;
vy4 = diff(y4,q1)*q_dot_1 + diff(y4,q2)*q_dot_2 + diff(y4,q3)*q_dot_3 + diff(y4,q4)*q_dot_4;

T4 = 0.5*m_4*[vx4 vy4]*[vx4; vy4;] + 0.5*I_4*((q_dot_1+q_dot_2+q_dot_3+q_dot_4)^2)

T = T1+T2+T3+T4

% Inertia Matrix
M11 = diff(diff(T, q_dot_1), q_dot_1);
M12 = diff(diff(T, q_dot_1), q_dot_2);
M13 = diff(diff(T, q_dot_1), q_dot_3);
M14 = diff(diff(T, q_dot_1), q_dot_4);

M21 = diff(diff(T, q_dot_2), q_dot_1);
M22 = diff(diff(T, q_dot_2), q_dot_2);
M23 = diff(diff(T, q_dot_2), q_dot_3);
M24 = diff(diff(T, q_dot_2), q_dot_4);

M31 = diff(diff(T, q_dot_3), q_dot_1);
M32 = diff(diff(T, q_dot_3), q_dot_2);
M33 = diff(diff(T, q_dot_3), q_dot_3);
M34 = diff(diff(T, q_dot_3), q_dot_4);

M41 = diff(diff(T, q_dot_4), q_dot_1);
M42 = diff(diff(T, q_dot_4), q_dot_2);
M43 = diff(diff(T, q_dot_4), q_dot_3);
M44 = diff(diff(T, q_dot_4), q_dot_4);

M = simplify([M11, M12, M13, M14; M21, M22, M23, M24; M31, M32, M33, M34; M41, M42, M43, M44])

% Coriolis and Centrifugal terms
q  = [q1; q2; q3; q4;];
dq = [q_dot_1; q_dot_2; q_dot_3; q_dot_4;];
C = sym(zeros(4,4));

for i = 1:4
    for j = 1:4
        tmp = 0;
        for k = 1:4
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
U2 = m_2*g0*(y2);
U3 = m_3*g0*(y3);
U4 = m_4*g0*(y4);

U = [U1+U2+U3+U4];
g_q = simplify([diff(U,q1); diff(U,q2); diff(U,q3); diff(U,q4)])