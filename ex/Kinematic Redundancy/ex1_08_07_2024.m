clear all
clc

%% === Variabili simboliche ===
syms alpha d d1 d2 d4 d6 l2 a1 a2 a3 a4 theta q1 q2 q3 q4 Dc3 a real
syms m1 m2 m3 I_c3xx I_c3yy I_c2yy I_c3zz real           % masse e momenti d'inerzia link
syms q1dot q2dot q3dot real                      % velocita' dei giunti
syms a1 a2 a3 a4 a5 real
% q1 => prismatico
% q2, q3 => rotazionali

%% number of joints 
N=3;

%% ============= Denavit-Hartenberg e calcolo cinematica diretta =================
DHTABLE = [ -pi/2   0     q1    0;     % giunto 1 prismatico
             pi/2   0     l2    q2;    % giunto 2 rotazionale
             0      Dc3   0     q3 ];  % giunto 3 rotazionale

TDH = [ cos(theta),               -sin(theta)*cos(alpha),    sin(theta)*sin(alpha),  a*cos(theta);
        sin(theta),                cos(theta)*cos(alpha),   -cos(theta)*sin(alpha),  a*sin(theta);
        0,                         sin(alpha),               cos(alpha),             d;
        0,                         0,                        0,                      1];

A = cell(1,N);
for i = 1:N
    alpha  = DHTABLE(i,1);
    a      = DHTABLE(i,2);
    d      = DHTABLE(i,3);
    theta  = DHTABLE(i,4);
    A{i}   = subs(TDH);
end

T = eye(4);
for i=1:N
    T = T*A{i};
    T = simplify(T);
end
T0N = T;  % Trasformazione finale dal frame 0 al frame N

%% ============= Posizione p_03 e sua Jacobiana analitica =================
A_0_1 = A{1};
A_0_2 = A{1}*A{2};
A_0_3 = simplify(A{1}*A{2}*A{3},'Steps',100);

p_03 = A_0_3(1:3,4);  % posizione del terzo link in base 0

px = p_03(1);
py = p_03(2);
pz = p_03(3);

dpx_dq1 = diff(px, q1); dpx_dq2 = diff(px, q2); dpx_dq3 = diff(px, q3);
dpy_dq1 = diff(py, q1); dpy_dq2 = diff(py, q2); dpy_dq3 = diff(py, q3);
dpz_dq1 = diff(pz, q1); dpz_dq2 = diff(pz, q2); dpz_dq3 = diff(pz, q3);

grad_px = [dpx_dq1; dpx_dq2; dpx_dq3];
grad_py = [dpy_dq1; dpy_dq2; dpy_dq3];
grad_pz = [dpz_dq1; dpz_dq2; dpz_dq3];

J_manual = [grad_px.'; grad_py.'; grad_pz.'];

disp('Jacobiana manuale = ')
disp(J_manual)
disp('Jacobiana manuale (semplificata) = ')
disp(simplify(J_manual,'Steps',100))

%% ============= Calcolo delle velocita' angolari con metodo geometrico ==========
R_0_1 = A_0_1(1:3,1:3);
R_0_2 = A_0_2(1:3,1:3);
R_0_3 = A_0_3(1:3,1:3);

z_0 = [0;0;1];
z_1_0 = R_0_1*z_0;  % asse giunto 2 in frame 0
z_2_0 = R_0_2*z_0;  % asse giunto 3 in frame 0

omega_1 = [0; 0; 0];                  % giunto 1 prismatico
omega_2 = omega_1 + q2dot*z_1_0;      % giunto 2 rotazionale
omega_3 = omega_2 + q3dot*z_2_0;      % giunto 3 rotazionale

disp(' ')
disp('=== Velocita angolari (Frame 0) ===')
disp('omega_1 = '), disp(omega_1)
disp('omega_2 = '), disp(omega_2)
disp('omega_3 = '), disp(omega_3)

omega_3_frame3 = R_0_3.' * omega_3;
disp('omega_3 in frame 3 = ')
disp(simplify(omega_3_frame3))

%% ============= Definizione Matrice di Inerzia M(q) =================
% Estrapolata dal tuo esempio, con i parametri (m1,m2,m3,I_c3xx, ecc.)
% e dipendente da (q2, q3). q1 e' prismatico, quindi non incide in M se 
% il link 1 e 2 non hanno contributi di rotazione. 
% NB: puoi adattare a seconda del modello completo.

M_q = [ m1+m2+m3,                   -m3*Dc3*cos(q2)*cos(q3),                      m3*Dc3*sin(q2)*sin(q3);
        -m3*Dc3*cos(q2)*cos(q3),    I_c3yy+I_c3xx+(I_c3yy - I_c3xx + m3*(Dc3^2))*(cos(q3)^2),     0;
         m3*Dc3*sin(q2)*sin(q3),    0,                                            I_c3zz+m3*(Dc3^2)];
     
M_q = simplify(M_q);

disp('Matrice di inerzia M(q) = ')
disp(M_q)

%% ============= Calcolo termini di Coriolis & centrifughi C(q,qd)*qd =============
% Formula standard (componente i-esima):
% c_i(q, qdot) = sum_{j,k} [1/2 * (dM_{i,j}/dq_k + dM_{i,k}/dq_j - dM_{j,k}/dq_i) ] * qdot_j * qdot_k

q  = [q1; q2; q3];      % giunti
qd = [q1dot; q2dot; q3dot];  % velocita'

% Inizializziamo il vettore c(q,qd) = [c1; c2; c3]
c_sym = sym(zeros(3,1));  

subsMap = [ m1+m2+m3,   I_c2yy+I_c3xx,  I_c3yy - I_c3xx + m3*(Dc3^2), ...
            m3*Dc3,   I_c3zz + m3*(Dc3^2) ];
varsMap = [ a1, a2, a3, a4, a5 ];

%% Sostituzione "dal basso verso l'alto":
% 1) Sostituisco i parametri fisici nei 5 coefficienti, in modo che a1..a5
%    diventino le "incognite"
M_in_terms_of_a = M_q;

% Nel codice seguente, useremo la funzione "issymmetric" per individuare
% corrispondenze dirette, ma è più semplice un passaggio manuale
M_in_terms_of_a = subs(M_in_terms_of_a, (m1+m2+m3), a1);             % a1
M_in_terms_of_a = subs(M_in_terms_of_a, (I_c2yy + I_c3xx), a2);       % a2
M_in_terms_of_a = subs(M_in_terms_of_a, (I_c3yy - I_c3xx + m3*Dc3^2), a3); % a3
M_in_terms_of_a = subs(M_in_terms_of_a, (m3*Dc3), a4);              % a4
M_in_terms_of_a = subs(M_in_terms_of_a, (I_c3zz + m3*Dc3^2), a5);    % a5

M_in_terms_of_a = simplify(M_in_terms_of_a);

disp(' ')
disp('=== Matrice di inerzia M(q) espressa con i coefficienti dinamici a1..a5 ===')
disp(M_in_terms_of_a)

for i=1:3
    expr = sym(0);
    for j=1:3
        for k=1:3
            expr = expr + 0.5 * ( diff(M_in_terms_of_a(i,j), q(k)) + diff(M_in_terms_of_a(i,k), q(j)) - diff(M_in_terms_of_a(j,k), q(i)) ) ...
                         * qd(j)*qd(k);
        end
    end
    c_sym(i) = simplify(expr);
end

disp('Vettore di Coriolis/centrifugo c(q,qd) = ')
disp(c_sym)
