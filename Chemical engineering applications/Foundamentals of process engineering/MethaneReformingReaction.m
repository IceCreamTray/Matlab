clc, clear all, close all, format long
                                                  %% REAZIONE DI STEAM REFORMING DEL GAS NATURALE

                                                                     % Legenda:
                                                               % CH4 + 2H2O = CO2 + 4H2
                                                               %  1      2      3    4 
                                                              
                                                              
                                                              

                                                          %% DEFINIZIONE DEI PARAMETRI

%% Condizioni iniziali
%inizializzazione delle variabili note
TCin=linspace(500,1200,36);    
disp(TCin)                               % vettore temperatura espresso in gradi celsius (da 500 a 1200 a passi di 20 gradi)
TKin=TCin+273.15;                        % vettore temperatura della miscela in entrata (espressa in gradi Kelvin)
P=10;                                    % Pressione cui avviene il processo (espressa in Bar)
Ppa= P*1e5;                              % Pressione cui avviene il processo (espressa in Pascal)

%% Costanti 
Rgas=8.314;                              % Costante universale dei gas ideali (unità di misura: Pascal*m^3/K/mol) 
Dhf=[-74.85 -241.83 -393.5 0];           % Entalpie di formazione standard diogni elemento, in ordine: CH4 H2O CO2 H2 (espresse in in kJ/mol)

%% Coefficienti dei calori specifici
%I coefficienti sono reperibili nelle tabelle. L'ordine per ciascun vettore è: CH4 H2O CO2 H2

a=[34.31*1e-3 33.46*1e-3 36.11*1e-3 28.84*1e-3];   
b=[5.469*1e-5 0.6880*1e-5 4.233*1e-5 0.00765*1e-5];
c=[0.3661*1e-8 0.7604*1e-8 -2.887*1e-8 0.3288*1e-8];
d=[-11.00*1e-12 -3.593*1e-12 7.464*1e-12 -0.8698*1e-12];

%Matrice dei coefficienti
C = [a; b; c; d];                        % Matrice 4x4: [4 coefficientix 4 composti]


%% Parametri cinetici e termodinamici
% Inizializzazione delle variabili cinetiche e termodinamiche note

k0=1.020e15;                              % fattore pre-esponeziale della costante cinetica espresso in kmol bar^0.5/g kgcat/h 
Eatt=243.9*1e3;                           % energia di attivazione della reazione espressa in J/mol
K0=6.194301277046940e11;                  % fattore pre-eesponenziale della costante di equilibrio espresso in bar^2
Dheq=192428.224210535;                    % Argomento dell'esponenziale della costante di equilibrio, espresso in J/mol


%% Scambi di calore

PowerQ=0;                                 %Potenza dovuta allo scambio termico, espressa in KW
Qdot=PowerQ*3600;                         %Calore scambiato dal reattore, espresso in kJ/h 


%% Geometria del reattore
%NOTA: il reattore viene considerato come volume di riferimento

D=0.4;                                    %Diametro, espresso in m
S=pi*D^2/4;                               %Sezione del reattore, espressa in m^2
v=0.1;                                    %Velocità del gas nel volume del reattore, espressa in m/s
Mcat=100;                                 %Massa del catalizzatore utilizzato, espressa in Kg
                                          
%% Valori molari per le specie in entrata
%NOTA: tutti i parametri molari in ingresso sono noti

Vin=S*v*3600;                             %portata volumetrica totale in ingresso,espressa in m^3/h
yin=[0.2 0.7 0.05 0.05];                  %frazioni molari dei gas in entrata: 20% di metano, 70% di acqua + prodotti residui.Sono grandezze adimensionali.


for i=1:1:length(TKin)


ntin=Ppa*Vin/Rgas/TKin(i);                   %portata molare totale in entrata, espressa in mol/h


%% Entalpie iniziali
%NOTA: viene chiamata una function
hin=calcolo_entalpia(1,TKin(i),C,Dhf);               %viene restituito un vettore contenente le entalpie di ciascuna specie, espresse in kJ/mol





                                                           %% METODO DI NEWTON-RAPHSON
                                                           
                                                           
%% Scelta della GUESS
y=[0.15 0.59 0.08];                          %Si ipotizzano le portate molari in uscita dei composti: CH4 H2o CO2. Sonograndezze adimensionali.
                                          %Si noti che la restante frazione molare è data da 1-y(19-y(2)-y(3)
nt=1.1;                                   %Si ipotizzano le moli totali in uscita, in forma adimensionale: nt=ntot_out/ntot_in
Tadim=0.86;                                %Si ipotizza la temperatura in uscita in forma adimensionale: Tadim=Tout/Tin
x0=[y(1),y(2),y(3), nt, Tadim]';          %Vettore colonna contenente la Guess (5 componenti). NOTA: Guess=incognite

%% Variabili necessarie al metodo
toll=1e-11;                               %Si definisce la tolleranza del metodo, errore minimo concesso
nmax = 100;                               %Si definisce il numero massimo di iterazioni possibili. 
                                          
                                          
                                          
%% Inizializzazione del numero di iterate e del vettore GUESS

n=0;                                     %numero iterate
x=x0;                                    %GUESS iniziale

%Entrata nel ciclo risolutivo
diff=toll+1;



                                                                %% Ciclo Risolutivo
while diff>=toll && n<nmax
    n=n+1;                                               %Conteggio del numero d'iterazioni
    
    %% Si riassegnano le variabili 
    %NOTA:sono tutte adimensionali
    y=x(1:3);                                           %Assegno le frazioni molari sono le prime 3 componenti del vettore GUESS
                                                        
    nt=x(4);                                            %Il numero di moli totali è la quarta componente del vettore GUESS                                           
    Tadim=x(5);                                         %la temperatura è la quinta componente del vettore GUESS
    y(4)=1-sum(y);                                      %La quarta frazione molare è uno meno la somma delle restanti. Pertanto è nota.
    
    %NOTA: ora il vettore y ha 4 componenti.
    
    %% Pressione parziale, velocità di reazione, costanti
    
    p=P*y;                                              %Vettore p contenente le pressioni parziali di ogni gas (costante P * vettore y)
    
    %Vengono definite le costanti
    k=k0*exp(-Eatt/(Rgas*TKin(i)*Tadim))*1e3;              %costante cinetica espressa in mol bar^0.5/kg cat/h 
    K=K0*exp(-Dheq/(Rgas*Tadim*TKin(i)));                  %costante di equilibrio espressa in bar^2 
    
    
    % velocità di reazione
    R=k*(p(1)*p(2)^2-p(4)^4*p(3)/K)/(p(4)^3.5);         %espressa in mol/Kg cat/h
    
    
    %% Calcolo delle entalpie e dei calori specifici
    
    %NOTA: vengono chiamate le function
    h=calcolo_entalpia (Tadim,TKin(i),C,Dhf);                     %entalpie in uscita espresse in KJ/mol
    cpi=calcolo_calore_specifico(Tadim,TKin(i),C);                         %calori specifici alla temperatura espressi in KJ/mol
     
    %% FUNZIONI
    
    %BILANCI DI MATERIA
    f(1)=yin(1) - nt*y(1) - R*Mcat/ntin;               % kmol/h - PER IL METANO
    f(2)=yin(2) - nt*y(2)-2*R*Mcat/ntin;               % kmol/h - PER L'ACQUA
    f(3)=yin(3) - nt*y(3)+R*Mcat/ntin;                 % kmol/h - PER L'ANIDRIDE CARBONICA
    f(4)=1-nt+2*R*Mcat/ntin;                           % kmol/h - BILANCIO TOTALE
    
    %BILANCIO DI ENERGIA
    f(5)=yin(1) + yin(2)*hin(2)/hin(1) + yin(3)*hin(3)/hin(1)+ ...
         yin(4)*hin(4)/hin(1)-y(1)*nt*h(1)/hin(1) - ...
         y(2)*nt*h(2)/hin(1)-y(3)*nt/hin(1)*h(3)- ...
         nt*y(4)/hin(1)*h(4)+Qdot/hin(1)/ntin;
     
    %Costruzione del vettore funzione (vettore colonna a 5 componenti)
    f=[f(1);f(2);f(3);f(4);f(5)];
    
    
    %% DERIVATE
    
    %%derivate di R
    dRy1=k*(((P^(-0.5)*y(2)^2)*((y(4)^(-3.5))+(y(1)*3.5*(y(4)^(-4.5)))))+(((P^1.5*y(3))/K)*0.5*(y(4)^(-0.5))));                    %derivata di R rispetto a y1        
    dRy2=k*(((P^(-0.5)*y(1))*((2*y(2)*(y(4)^(-3.5)))+(y(2)^2*3.5*(y(4)^(-4.5)))))+(((P^1.5*y(3))/K)*(0.5*(y(4)^(-0.5)))));         %derivata di R rispetto a y2
    dRy3=k*((P^(-0.5)*y(1)*y(2)^2*3.5*(y(4)^(-4.5)))-((P^1.5/K)*((y(4)^0.5)-(y(3)*0.5*(y(4)^(-0.5))))));                           %derivata di R rispetto a y3
    
    dRTadim=1e3*(k0*(Rgas*TKin(i)*Tadim^2*P^0.5)^(-1))*(y(4)^0.5)*exp(-Eatt/(Rgas*Tadim*TKin(i)))...
        *((y(1)*y(2)^2*Eatt/(y(4)^(4)))-((P^2*y(3)/K0)*(Eatt-Dheq)*exp(Dheq*(Rgas*TKin(i)*Tadim)^(-1))));                             %derivata di R rispetto a Tadim
    
 



    %derivate delle entalpie di ciascuna specie rispetto alla temperatura
    dh1dt=cpi(1)*TKin(i);
    dh2dt=cpi(2)*TKin(i);
    dh3dt=cpi(3)*TKin(i);
    dh4dt=cpi(4)*TKin(i);
    
    %Costruzione delle derivate delle funzioni rispetto ad ogni componente
    
    %derivate parziali di f1
    df11=-nt-(Mcat/ntin)*dRy1;                                                 
    df12=-(Mcat/ntin)*dRy2;                                                    
    df13=-(Mcat/ntin)*dRy3;
    df14=-y(1);
    df15=-(Mcat/ntin)*dRTadim;

    %derivate parziali di f2
    df21=-(2*(Mcat/ntin))*dRy1;
    df22=-nt-(2*(Mcat/ntin))*dRy2;
    df23=-(2*(Mcat/ntin))*dRy3;
    df24=-y(2);
    df25=-(2*(Mcat/ntin))*dRTadim;

    
    %derivate parziali di f3
    df31=(Mcat/ntin)*dRy1;
    df32=(Mcat/ntin)*dRy2;
    df33=-nt+(Mcat/ntin)*dRy3;
    df34=-y(3);
    df35=(Mcat/ntin)*dRTadim;

    %derivate parziali di f4
    df41=(2*(Mcat/ntin))*dRy1;
    df42=(2*(Mcat/ntin))*dRy2;
    df43=(2*(Mcat/ntin))*dRy3;
    df44=-1;
    df45=(2*(Mcat/ntin))*dRTadim;

    
    %derivate parziali di f5
    df51=-nt*(h(1)/hin(1))+nt*(h(4)/hin(1));
    df52=-nt*(h(2)/hin(1))+nt*(h(4)/hin(1));
    df53=-nt*(h(3)/hin(1))+nt*(h(4)/hin(1));
    df54=-(y(1)*(h(1)/hin(1)))-(y(2)*(h(2)/hin(1)))-(y(3)*(h(3)/hin(1)))-((y(4))*(h(4)/hin(1)));
    df55=-((y(1)*(nt/hin(1)))*dh1dt)-((y(2)*(nt/hin(1)))*dh2dt)-((y(3)*(nt/hin(1)))*dh3dt)-((y(4)*(nt/hin(1)))*dh4dt); 
    
    %% Risoluzione

    %Costruzione dellamatrice Jacobiana (5x5)
     df=[df11 df12 df13 df14 df15
         df21 df22 df23 df24 df25
         df31 df32 df33 df34 df35
         df41 df42 df43 df44 df45
         df51 df52 df53 df54 df55];
    
    %controllo sul determinante
    if det(df)== 0
         disp ('errore: determinante nullo'); 
    end
    
    %Aggiornamento delle variabili
    diff=-df\f;                                       
    x=x+diff;
    
    %stima dell'errore come norma a infinito
    diff=max(abs(diff));
    
    %controllo sul numero delle iterate
     if n==nmax
       disp ('errore: numero massimo di iterate raggiunto') 
       return
    end
    
end



%conversione del ch4: 
X(i)=(yin(1)-x(1))/yin(1);

% y3(i)=x(3);                frazioni molari prodotti
% y4(i)=(1-x(3)-x(2)-x(1));


disp('frazione molare di metano:');
disp(x(1));
disp('frazione molare di acqua');
disp(x(2));
disp('frazione molare di anidride carbonica');
disp(x(3));
disp('frazione molare di idrogeno');
disp((1-x(1)-x(2)-x(3)));
disp('moli in uscita (adim)');
disp(x(4));
disp('Temperatura in uscita (adim)');
disp(x(5));
disp('numero di iterate');
disp(n);

end



plot(TCin,X);  %plotta il grafico conversione temperatura
% plot(TCin,y3)%plotta il grafico frazioni molari-temperatura
% hold on;
% plot(TCin,y4)
















