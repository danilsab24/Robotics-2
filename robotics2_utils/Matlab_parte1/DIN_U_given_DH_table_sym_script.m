%numero di joint
n=3;

syms rcx1 rcx2 rcx3 rcy1 rcy2 rcy3 real;


q = sym('q', [n, 1]);
for i=1:n
    nomeVariabile = ['q' num2str(i)];
    assignin('base', nomeVariabile,q(i,1));
end

L = sym('l', [n, 1]);
for i=1:n
    nomeVariabile = ['l' num2str(i)];
    assignin('base', nomeVariabile,L(i,1));
end

M = sym('m', [n, 1]);
for i=1:n
    nomeVariabile = ['m' num2str(i)];
    assignin('base', nomeVariabile,M(i,1));
end

rcx = sym('rcx', [1,n]);
for i=1:n
    nomeVariabile = ['rcx' num2str(i)];
    assignin('base', nomeVariabile,rcx(1,i));
end
rcy = sym('rcy', [1,n]);
for i=1:n
    nomeVariabile = ['rcy' num2str(i)];
    assignin('base', nomeVariabile,rcy(1,i));
end
rcz = sym('rcz', [1,n]);
for i=1:n
    nomeVariabile = ['rcz' num2str(i)];
    assignin('base', nomeVariabile,rcz(1,i));
end

rcz=subs(rcz,rcz,[0,0,0])  % da cancellare se si vogliono fare calcoli con rcz

uni=ones(1,n)
com=[rcx;rcy;rcz;uni]


%DH_table
DH_table=[0 l1 0 q1; 0 l2 0 q2; 0 l3 0 q3 ]

%calcolo delle Denavit-Hartemberg Matrices
syms a alpha tetha d real;
DH_param=[alpha a d tetha];

DH_MAT=[cos(tetha) -cos(alpha)*sin(tetha) sin(alpha)*sin(tetha)  a*cos(tetha);
        sin(tetha)  cos(alpha)*cos(tetha) -sin(alpha)*cos(tetha) a*sin(tetha);
        0           sin(alpha)              cos(alpha)                d      ;
        0               0                       0                     1      ];

Array_A = sym('A', [4, 4, n]);
for i=1:n
    DH_MAT_to_sub=subs(DH_MAT,DH_param,DH_table(i,:));
    Array_A(:,:,i)=DH_MAT_to_sub;
end

Array_A;


u = sym('u', [n, 1]);
for i=1:n
    nomeVariabile = ['u' num2str(i)];
    assignin('base', nomeVariabile,u(i,1));
end

syms g G real
G=[0, g, 0 ,0]

for i=1:n
    partial=eye(4);
    for j=1:i
        partial=simplify(partial*Array_A(:,:,j));
    end
    u(i,1)=simplify(-M(i,1)*G*partial*com(:,i));
end
u

U=u(1,1);
for i=2:n
    U=U+u(i,1);
end
U


G_1=simplify(diff(U,q1));
G_2=simplify(diff(U,q2));
G_3=simplify(diff(U,q3));
G=[G_1;G_2;G_3]


eq= G_3==0
solve(eq,[l3,rcx3,rcy3])

simplify(subs(G_2,[l3,rcy3],[-rcx3,0]))
%simplify(subs(G_1),[rcx2,rcy2],[-((m3+m2)/m2)*l2,0])
