clc

% pxdes=2;
% pydes=1;
% pzdes=


% syms a b c;  % Definisci le variabili simboliche
%         eq1 = 1 == a+b+c
%         eq2 = 0 == 3*a+2*b+c
%         eq3 = 0 == 6*a+2*b
% 
% solutions = solve([eq1, eq2, eq3], [a, b, c])


syms a b ;  % Definisci le variabili simboliche
        eq1 = 1 == a+b
        eq2 = 0 == -3/2*a-1*b


solutions = solve([eq1, eq2], [a, b])

% syms theta1 theta2 theta3;  % Definisci le variabili simboliche
%         %eq1 = pxdes == K*cos(q1)-q2*sin(q1)+D*cos(q1+q3)
%         %eq2 = pydes == K*sin(q1)+q2*cos(q1)+D*sin(q1+q3)
%         %eq3 = phides == q1+q3
%         
%         eq1= pxdes == a1*cos(theta1) + d2*sin(theta1) + a2*cos(theta1)*cos(theta2)
%         eq2= pydes == a1*sin(theta1) - d2*cos(theta1) + a2*cos(theta2)*sin(theta1)
%         eq3= pzdes ==        a2*sin(theta2)
% 
% 
% 
% solutions = solve([eq1, eq2, eq3], [theta1, theta2, theta3]);





% disp('Le soluzioni per q1 sono:');
% disp(double(solutions.q1))
% 
% disp('Le soluzioni per q2 sono:');
% disp(double(solutions.q2))
% 
% disp('Le soluzioni per q3 sono:');
% disp(double(solutions.q3))
