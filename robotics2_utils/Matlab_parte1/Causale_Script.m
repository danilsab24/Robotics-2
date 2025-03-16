% n=4; 
% q = sym('q', [n, 1]);
% for i=1:n
%     nomeVariabile = ['q' num2str(i)];
%     assignin('base', nomeVariabile,q(i,1));
% end
% 
% L = sym('l', [n, 1]);
% for i=1:n
%     nomeVariabile = ['l' num2str(i)];
%     assignin('base', nomeVariabile,L(i,1));
% end
% 
% %inertia matrix for each link
% I = sym('Ic', [3, 3, n]);
% for i=1:n
%     I(1,1,i)=sym(['Ic' num2str(i) 'xx']);I(1,2,i)=sym(['Ic' num2str(i) 'xy']); I(1,3,i)=sym(['Ic' num2str(i) 'xz']);
%     I(2,1,i)=sym(['Ic' num2str(i) 'yx']);I(2,2,i)=sym(['Ic' num2str(i) 'yy']); I(2,3,i)=sym(['Ic' num2str(i) 'yz']);
%     I(3,1,i)=sym(['Ic' num2str(i) 'zx']);I(3,2,i)=sym(['Ic' num2str(i) 'zy']); I(3,3,i)=sym(['Ic' num2str(i) 'zz']);
% end
% 
% M = sym('m', [n, 1]);
% for i=1:n
%     nomeVariabile = ['m' num2str(i)];
%     assignin('base', nomeVariabile,M(i,1));
% end
% 
% pcx = sym('pcx', [1,n]);
% for i=1:n
%     nomeVariabile = ['pcx' num2str(i)];
%     assignin('base', nomeVariabile,pcx(1,i));
% end
% pcy = sym('pcy', [1,n]);
% for i=1:n
%     nomeVariabile = ['pcy' num2str(i)];
%     assignin('base', nomeVariabile,pcy(1,i));
% end
% pcz = sym('pcz', [1,n]);
% for i=1:n
%     nomeVariabile = ['pcz' num2str(i)];
%     assignin('base', nomeVariabile,pcz(1,i));
% end
% 
% pcz=subs(pcz,pcz,zeros(1,n));  % da cancellare se si vogliono fare calcoli con pcz
% 
% uni=ones(1,n);
% com=[pcx;pcy;pcz;uni]
% 
% 
% D = sym('dc', [n, 1]);
% for i=1:n
%     nomeVariabile = ['dc' num2str(i)];
%     assignin('base', nomeVariabile,D(i,1));
% end
% 
% 
% %assegnazioni rispetto i lframe di base
% %assegnamo la prima posizione del centro di massa:
% pcx1= dc1*(cos(q1));
% pcy1=dc1*(sin(q1));
% pcz1=0; %era zero e rimane zero 
% 
% %assegnamo la prima posizione del centro di massa:
% pcx2= l1*(cos(q1))+dc2*cos(q2);
% pcy2=l1*(sin(q1))-dc2*sin(q2);
% pcz2=0; %era zero e rimane zero 
% 
% %assegnamo la prima posizione del centro di massa:
% pcx3= l1*(cos(q1))+l2*cos(q2)+dc3*cos(q3);
% pcy3=l1*(sin(q1))-l2*sin(q2)+dc3*sin(q3);
% pcz3=0; %era zero e rimane zero 
% 
% %assegnamo la prima posizione del centro di massa:
% pcx4=l1*(cos(q1))+l2*cos(q2)+l3*cos(q3)+dc4*cos(q4);
% pcy4=l1*(sin(q1))-l2*sin(q2)+l3*sin(q3)+dc4*sin(q4);
% pcz4=0; %era zero e rimane zero 
% 
% com(1,1)




J=[1 1 3/2 1 ; sqrt(3)+1 sqrt(3) sqrt(3)/2 0]