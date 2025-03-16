clear all
clc

%% Define symbolic variables

syms alpha d a theta

%% number of joints of Robot

disp('Inserisci il numero --INTERO-- di Joint: ');
N=chiediNumero();

%% Insert DH table of parameters of RObot Con possibilità di correggere la riga3

DHTABLE= cell(N,4);
insertRow=cell(1,4);
typeOfJoints=char(N,1);
tipoJoint=' ';
q=sym('q',[1,N]);

for i=1:N
    while true
        disp([10 10 'Inserisci la riga n.' num2str(i) ' della tabella DH.']);
        tipoJoint=chiediTipoJoint(i);
        for j=1:4
            result=richiediparametro(i,j);
            if result<=-1000
                if j==1  param='alpha';  end
                if j==2  param='a';  end
                if j==3  param='d';  end
                if j==4  param='theta';  end
                nomeVariabile = [param num2str(i)];
                resultSym = sym(nomeVariabile);
                % Crea la variabile nello spazio di lavoro principale
                assignin('base', nomeVariabile, resultSym);
                insertRow{1,j}=resultSym;
                if result==-5000
                     q(1,i)=resultSym;
                end
            else
                insertRow{1,j}=result;
            end
        end
        
        disp([10 'Hai inserito questa riga n. ' num2str(i) ':' ]);
        disp(insertRow);
        disp("E' corretta? Attenzione anche al tipo di Joint!!!");
        
        RigaCorretta=chiediSiNo();

        if(RigaCorretta==true) 
            for k=1:4
                DHTABLE{i,k}=insertRow{1,k};
                typeOfJoints(i,1)=tipoJoint;
            end
            break;          
        end
    end
end
disp([10 'D-H table given in input:'])
disp(DHTABLE);



%For Scara
%DHTABLE = [  0   sym('a1') sym('d1') sym('q1');
%             0   sym('a2')    0      sym('q2');
%             0     0       sym('d3')    0;
%             pi    0       sym('d4') sym('q4')];
% 
% DHTABLE = [  0   sym('a1') sym('d1') sym('q1');
%             pi/2   0 sym('d2')  sym('q2');
%             0     sym('a3')   0   sym('q3')   ;
%             0     sym('a4') 0 sym('q4')];



         
%% Build the general Denavit-Hartenberg trasformation matrix

TDH = [ cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta);
        sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);
          0             sin(alpha)             cos(alpha)            d;
          0               0                      0                   1];


Rot_Mat=[ cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha); 
        sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) ;
          0             sin(alpha)             cos(alpha)        ];


%% Build transformation matrices for each link
% First, we create an empty cell array

A=cell(1,N);
R=cell(1,N);
matRot=sym('matRot',[3,3,N]);
zetas=cell(1,N);


% For every row in 'DHTABLE' we substitute the right value inside
% the general DH matrix

for i = 1:N
    alpha = DHTABLE(i,1);
    a = DHTABLE(i,2);
    d = DHTABLE(i,3);
    theta = DHTABLE(i,4);
    A{i} = subs(TDH);
end

for i = 1:N
    alpha = DHTABLE(i,1);
    a = DHTABLE(i,2);
    d = DHTABLE(i,3);
    theta = DHTABLE(i,4);
    R{i} = subs(Rot_Mat);
end


%% Direct kinematics
disp([10 10 'Direct kinematics of robot in symbolic form (simplifications may need some time)'])
disp(['Number of joints N=' num2str(N) 10 10])


for i = 1:N
    disp(['matrice  ' num2str(i-1) 'A' num2str(i) '= ']);
    disp(A{i});
    %det(A{i})
end

% for i = 1:N
%     disp(['matrice  ' num2str(i-1) 'R' num2str(i) '= ']);
%     disp(R{i});
%     %det(R{i})
% end

matRot(:,:,1)=R{1};
for i=2:N
    matRot(:,:,i)=simplify(matRot(:,:,i-1)*R{i});
end


z0=[0;0;1];
k=0;
for i= 1:N
    r=eye(3);
    for j=1:k  
        r=r*R{j};
        r = simplify(r);
    end
    zetas{i}=r*z0;
    k=k+1;
end


% Note: 'simplify' may need some time
T = eye(4)
for i=1:N 
    T = T*A{i}
    T = simplify(T);
end
%det(T)


% output TN matrix
T0N = T
% output ON position
p = T(1:3,4)
disp('note that the same result can be computed in a efficient way by: ');
disp('this formula 0A1*(1A2*(2A3*(..*(jAn * [0; 0; 0; 1]) ) ) )')
disp('dove [0; 0; 0; 1] è la posizione dell originie dell end effector');
disp('nella terna di riferimento dell end effector');
disp('this way is more efficient when we have 3+ joint');
% output xN axis
n=T(1:3,1)
% output yN axis
s=T(1:3,2)
% output zN axis
a=T(1:3,3)

J = simplify(jacobian(p, q),Steps=100);
J_Geo=Build_Jacobian(J,zetas,typeOfJoints,N);
disp(J_Geo)
J_geo_nuovo=cambioRiferimentoJacobiano(J_Geo,matRot)
sol=checkSingularities(J_geo_nuovo,q)
%disp(sol)


% disp('Le soluzioni per var1 sono:');
% disp(double(sol.theta1))
% 
% disp('Le soluzioni per var2 sono:');
% disp(double(sol.theta2))
% 
% disp('Le soluzioni per var3 sono:');
% disp(double(sol.d3))

% disp('Le soluzioni per var4 sono:');
% disp(double(sol.theta4))



%sostituzione dei valori delle variabili dentro J_Geo
%J_Geo_Val=subs(J_Geo,q,[0,0,0,0])






















function int=chiediNumero()
    %chiede funche non metti un numero corretto, con exit termini il codice
    while true
        % Chiedi all'utente di inserire una stringa
        userInputString = input('Inserisci un numero, con "exit" termina il codice: ', 's');
            
        if strcmp(userInputString, 'exit')
             error('Terminazione anticipata del codice.');
        end

        %Vediamo se la stringa è un espressione matematica
        try
            risultato = eval(userInputString);
            disp(['Hai inserito il numero: ' num2str(risultato)]);
            int=risultato;
            return;
        catch
            disp('L input non è un espressione matematica');
        end
        %Vediamo se la stringa inserita è un numero 
        numeroInserito = str2double(userInputString);
        if ~isnan(numeroInserito)
            disp(['Hai inserito il numero: ' num2str(numeroInserito)]);
            int=numeroInserito;
            return;  % Restituisci il numero se la conversione ha avuto successo
        else
            disp("L' input non è un numero valido. Riprova.");
        end
    end
end








function J_nuovo=cambioRiferimentoJacobiano(J,matRot)
    disp(['Vuoi cambiare il riferimento in cui è espresso lo Jacobiano? ' 10]);
    if(chiediSiNo)
        while true
            userInputString = input('In che riferimento vuoi esprimere lo jacobiano ? :', 's');
            if strcmp(userInputString, 'exit')
                 error('Terminazione anticipata del codice.');
            end
            numeroInserito = str2double(userInputString);
            if ~isnan(numeroInserito)
                disp(['Hai inserito il numero: ' num2str(numeroInserito)]);
                int=numeroInserito;
                break; 
            else
                disp("L' input non è un numero valido. Riprova.");
            end
        end
        
        M=[transpose(matRot(:,:,int)),zeros(3,3);zeros(3,3),transpose(matRot(:,:,int))];
        J_nuovo=M*J;
        J_nuovo=simplify(J_nuovo);
    else
        J_nuovo=J;
    end
end



function sol=checkSingularities(J,q)
    dim=size(J);
    if(dim(1,1)==dim(1,2))
        eqn=det(J)==0;
        sol = solve(eqn,q);
    elseif(dim(1,1)>dim(1,2))
        disp('det non semplificato');
        det(transpose(J)*J)
        daRis=simplify(det(transpose(J)*J),'Steps',50);
        eqn= daRis==0
        sol=solve(eqn,q)
        disp(sol);
    else
        disp('det non semplificato');
        det(J*transpose(J))
        daRis=simplify(det(J*transpose(J)),'Steps',50);
        eqn= daRis ==0
        sol=solve(eqn,q)
        disp(sol);
    end
end


function J_Geo=Build_Jacobian(J,zetas,typeOfJoints,N)
    J_Geo=sym('J_Geo',[6,N]);
    z=[0;0;0];
    for i=1:3
        for j=1:N
            J_Geo(i,j)=simplify(J(i,j));
        end
    end
    for i=1:N
        if typeOfJoints(i)=='R'
            z=zetas{i};
            for j=4:6
                J_Geo(j,i)=z(j-3);
            end
       else
            for j=4:6
                J_Geo(j,i)=0;
            end
        end
    end
end



function tipoJoint=chiediTipoJoint(indice)
    while(true)
        userInputString = input('Inserisci il tipo di Joint: ' , 's');
        if (strcmp(userInputString, 'R')|| strcmp(userInputString, 'r')  )
            tipoJoint='R';
            disp(['il giunto ' num2str(indice) ' è Rotazionale']);
            return;
        elseif (strcmp(userInputString, 'P')|| strcmp(userInputString, 'p')  )  
                tipoJoint='P';
                disp(['il giunto ' num2str(indice) ' è Prismatico']);
                return;
        elseif strcmp(userInputString, 'exit')
                error('Terminazione anticipata del codice.');
        else
            disp('Inserire una tipologia di giunti valida');
        end
    end
end


function parametroInserito = richiediparametro(indice,tipoParametro)
    while true
        % Chiedi all'utente di inserire una stringa
        userInputString = input(['elem(' num2str(indice) ',' num2str(tipoParametro) ']. ' 'Inserisci un numero o "var" le variabili: ' ], 's');
        if strcmp(userInputString, 'var')  
            if tipoParametro==1  param='alpha';  end
            if tipoParametro==2  param='a';  end
            if tipoParametro==3  param='d';  end
            if tipoParametro==4  param='theta';  end
            disp(['Hai inserito la variabile: ' param  num2str(indice)]);
            parametroInserito=-5000;
            return;

        elseif strcmp(userInputString, 'pa')  
            if tipoParametro==1  param='alpha';  end
            if tipoParametro==2  param='a';  end
            if tipoParametro==3  param='d';  end
            if tipoParametro==4  param='theta';  end
            disp(['Hai inserito il parametro: ' param  num2str(indice)]);
            parametroInserito=-1000;
            return;
        else
            if strcmp(userInputString, 'exit')
                error('Terminazione anticipata del codice.');
            end

            %Vediamo se la stringa è un espressione matematica
            try
                risultato = eval(userInputString);
                disp(['Hai inserito il numero: ' num2str(risultato)]);
                parametroInserito=risultato;
                return;
            catch
                disp('L input non è un espressione matematica');
            end
            %Vediamo se la stringa inserita è un numero 
            numeroInserito = str2double(userInputString);
            if ~isnan(numeroInserito)
                disp(['Hai inserito il numero: ' num2str(numeroInserito)]);
                parametroInserito=numeroInserito;
                return;  % Restituisci il numero se la conversione ha avuto successo
            else
                disp("L' input non è un numero valido. Riprova.");
            end
        end
    end
end

function bool=chiediSiNo()
    bool=false;
    %chiede funche non metti un numero corretto, con exit termini il codice
    while true
        % Chiedi all'utente di inserire una stringa
        userInputString = input('Rispondi si/no : ', 's');
        if strcmp(userInputString, 'si')bool=true;return; end
        if strcmp(userInputString, 's')bool=true;return; end
        if strcmp(userInputString, 'SI')bool=true;return; end
        if strcmp(userInputString, 'Si')bool=true;return; end
        if strcmp(userInputString, 'sI')bool=true;return; end
        if strcmp(userInputString, 'yes')bool=true;return; end
        if strcmp(userInputString, 'y')bool=true;return; end
        if strcmp(userInputString, 'YES')bool=true;return; end
        if strcmp(userInputString, 'Yes')bool=true;return; end
        if strcmp(userInputString, 'YeS')bool=true;return; end
        if strcmp(userInputString, 'YEs')bool=true;return; end
        if strcmp(userInputString, 'NO')bool=false;return; end
        if strcmp(userInputString, 'no')bool=false;return; end
        if strcmp(userInputString, 'n')bool=false;return; end
        if strcmp(userInputString, 'nO')bool=false;return; end
        if strcmp(userInputString, 'No')bool=false;return; end
        if(strcmp(userInputString, 'exit') ) error('Terminazione anticipata del codice.'); end
    end
end
