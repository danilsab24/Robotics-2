clear all
clc

%% Robotics 2 Midterm
% 24 April 2024 - Ex #1
%
% Projected Gradient (PG) method (with obstacle avoidance)
% and Task Priority (TP) method (with two variants)
% for a 3R planar robot

disp('*** Projected Gradient and Task Priority (TP) methods for a 3R planar robot ***');
disp(' ');

syms q1 q2 q3 real
q = [q1; q2; q3];

disp('Kinematics of the 3R planar robot (links of unitary length)');
pe_simb = [cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
           sin(q1)+sin(q1+q2)+sin(q1+q2+q3)];
pm_simb = [cos(q1)+cos(q1+q2);
           sin(q1)+sin(q1+q2)];

Je_simb = jacobian(pe_simb,q);
Jm_simb = jacobian(pm_simb,q);

disp('At the given configuration (as in the text figure)');
q0 = [0; pi/2; -pi/2];

disp('Numerical kinematics');
pe = subs(pe_simb,q,q0);
pm = subs(pm_simb,q,q0);
Je = subs(Je_simb,q,q0);
Jm = subs(Jm_simb,q,q0);

disp('Clearance');
C = [0; 2];
r = 0.5;
H = norm(pm - C) - r

disp('Desired end-effector velocity');
ve = [0; 1];

disp('Projected Gradient (PG) method');
nablaH = 0.5 * Jm' * (pm - C) / norm(pm - C);
nablaH = eval(nablaH);
alfa1 = 1; % Funziona anche con valori diversi da 1
Je_pinv = pinv(Je);

% dq_PG = alfa1*nablaH + Je_pinv*(ve - alfa1*Je*nablaH);
dq_PG = alfa1*nablaH + Je_pinv*(ve - alfa1*(Je*nablaH));
dq_PG = eval(dq_PG);

disp('Check PG solution');
ve_check = eval(Je * dq_PG)
vm = eval(Jm * dq_PG)

% TP con due casi: A & B
disp('Task Priority (TP) - case A');
Pe = eye(3) - Je_pinv * Je;
vm_A = vm;
dq_TP_A = Je_pinv * ve + pinv(Jm * Pe) * (vm_A - Jm * (Je_pinv * ve));
dq_TP_A = eval(dq_TP_A);

disp('Check TP A solution');
ve_check = eval(Je * dq_TP_A)
vm_check = eval(Jm * dq_TP_A)

disp('Task Priority (TP) - case B');
alfa2 = 1; % Con alfa2=0.9 si minimizza l'errore vm - Jm*dq_TP_B
vm_B = alfa2 * (1 - (r / norm(pm - C))) * (pm - C);
vm_B = eval(vm_B);

dq_TP_B = Je_pinv * ve + pinv(Jm * Pe) * (vm_B - Jm * (Je_pinv * ve));
dq_TP_B = eval(dq_TP_B);

disp('Check TP B solution');
ve_check = eval(Je * dq_TP_B)
vm_check = eval(Jm * dq_TP_B)
vm_check_error = eval(vm_B - Jm * dq_TP_B)
