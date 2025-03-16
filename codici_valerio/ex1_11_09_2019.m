clear all
clc

syms alpha d d1 d2 d4 d6 a l1 l2 l3 a4 theta q1 q2 q3 q4 a2 a3 

%% number of joints 
N=3;


%% PAY ATTENTION TO THE POSITION OF
%% a and d: a is the second column
%% d the third!

DHTABLE = [  pi/2   0    d1    q1;
             0      a2   0     q2;
             0      a3   0     q3];

         
TDH = [ cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta);
        sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);
          0             sin(alpha)             cos(alpha)            d;
          0               0                      0                   1];

A = cell(1,N);

for i = 1:N 
    alpha = DHTABLE(i,1);
    a = DHTABLE(i,2);
    d = DHTABLE(i,3);
    theta = DHTABLE(i,4);
    A{i} = subs(TDH);
    disp(i)
    disp(A{i})
end


T = eye(4);

for i=1:N 
    T = T*A{i};
    T = simplify(T);
end

T0N = T;

p = T(1:3,4);

n = T(1:3,1);

s = T(1:3,2);

a = T(1:3,3);

A_0_1 = A{1};

A_0_2 = A{1} * A{2};

A_0_3 = simplify(A{1} * A{2} * A{3},steps=100)

p_01 = A_0_1(1:3, end);
p_02 = A_0_2(1:3, end);
p_03 = A_0_3(1:3, end);

p_0_E = p_03;
p_1_E = p_03 - p_01;
p_2_E = p_03 - p_02;
p_3_E = p_03 - p_03

R_0_1 = A_0_1(1:3, 1:3);
R_0_2 = A_0_2(1:3, 1:3);
R_0_3 = A_0_3(1:3, 1:3)

z_0 = [0;
       0;
       1];

z_1 = simplify(R_0_1*z_0, Steps=100);
z_2 = simplify(R_0_2*z_0, Steps=100);

p_z_0 = simplify(cross(z_0, p_0_E), steps = 100);
p_z_1 = simplify(cross(z_1, p_1_E), steps = 100);
p_z_2 = simplify(cross(z_2, p_2_E), steps = 100);


%% TO MODIFY BASED ON THE TYPE OF JOINTS
J_L_A = [p_z_0, p_z_1, p_z_2;
         z_0, z_1, z_2]


rotation_mid_frame_1 = [R_0_1, [0, 0, 0; 
                                  0, 0, 0; 
                                  0, 0, 0;];
                        [0, 0, 0; 
                         0, 0, 0; 
                         0, 0, 0;], R_0_1];


rotation_mid_frame_2 = [R_0_2, [0, 0, 0; 
                                  0, 0, 0; 
                                  0, 0, 0;];
                        [0, 0, 0; 
                         0, 0, 0; 
                         0, 0, 0;], R_0_2];



rotation_mid_frame_3 = [R_0_3, [0, 0, 0; 
                                  0, 0, 0; 
                                  0, 0, 0;];
                        [0, 0, 0; 
                         0, 0, 0; 
                         0, 0, 0;], R_0_3];

J_l_subs = double(subs(J_L_A(1:3, 1:3), [a2,a3,d1,q1,q2,q3],[4,3,5,pi/2,pi/4,pi/2]))

det_J_l_subs = det(J_l_subs)

syms  q1_dot q2_dot q3_dot q1(t) q2(t) q3(t) 

J_l = [ -sin(q1(t)) * (a3*cos(q2(t) + q3(t)) + a2*cos(q2(t))),  -cos(q1(t)) * (a3*sin(q2(t) + q3(t)) + a2*sin(q2(t))),  -a3*sin(q2(t) + q3(t))*cos(q1(t));
       cos(q1(t)) * (a3*cos(q2(t) + q3(t)) + a2*cos(q2(t))),  -sin(q1(t)) * (a3*sin(q2(t) + q3(t)) + a2*sin(q2(t))),  -a3*sin(q2(t) + q3(t))*sin(q1(t));
       0,                                                     a3*cos(q2(t) + q3(t)) + a2*cos(q2(t)),                 a3*cos(q2(t) + q3(t))];

J_l_dot = diff(J_l, t)

n_q_q_do_mat = J_l_dot*[q1_dot; q2_dot; q3_dot]

syms q1_dot q2_dot q3_dot q1 q2 q3 a2 a3 d1
J_l_dot_mat = [ sin(q1)*(a3*sin(q2 + q3)*(q2_dot + q3_dot) + a2*sin(q2)*q2_dot) - cos(q1)*(a2*cos(q2) + a3*cos(q2 + q3))*q1_dot,   sin(q1)*q1_dot*(a2*sin(q2) + a3*sin(q2 + q3)) - cos(q1)*(a3*cos(q2 + q3)*(q2_dot + q3_dot) + a2*cos(q2)*q2_dot),   a3*sin(q1)*sin(q2 + q3)*q1_dot - a3*cos(q1)*cos(q2 + q3)*(q2_dot + q3_dot);
                 - cos(q1)*(a3*sin(q2 + q3)*(q2_dot + q3_dot) + a2*sin(q2)*q2_dot) - sin(q1)*(a2*cos(q2) + a3*cos(q2 + q3))*q1_dot, - sin(q1)*(a3*cos(q2 + q3)*(q2_dot + q3_dot) + a2*cos(q2)*q2_dot) - cos(q1)*(a2*sin(q2) + a3*sin(q2 + q3))*q1_dot, - a3*cos(q1)*sin(q2 + q3)*q1_dot - a3*sin(q1)*cos(q2 + q3)*(q2_dot + q3_dot);
                 0, - a3*sin(q2 + q3)*(q2_dot + q3_dot) - a2*sin(q2)*q2_dot, -a3*sin(q2 + q3)*(q2_dot + q3_dot)];

J_l_dot_mat_subs = double(subs(J_l_dot_mat,[q1,q2,q3,q1_dot,q2_dot,q3_dot,a2,a3],[pi/2,pi/4,pi/2,1,2,-2,4,3]))

n_q_q_do_mat_subs = J_l_dot_mat_subs*[1;2;-2]

q_dot_dot = inv(J_l_subs)*n_q_q_do_mat_subs