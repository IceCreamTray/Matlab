%% CALCOLO DELLA TEMPERATURA DELL'ACQUA DI RINVERDIMENTO IN USCITA DAL CONDENSATORE


%% DEFINIZIONE DEI DATI
Q=406953; %potenza termica scambiata (J/s)
MH2O=20; %portata di acqua (kg/s)
MWH2O=0.018; %peso molecolare acqua (kg)
Tin=16; %temperatura dell'acqua in entrata (Celsius)
Tout=30; %guess (celsius)
diff=1;
%% DEFINIZIONE COSTANTI
c1=276370; 
c2=-2090.1; 
c3=8.125; 
c4=-0.014116; 
c5=9.3701*1e-06;


%% DEFINIZIONE VETTORI
x=[];
y=[];
z=[];

%% ENTRATA NEL LOOP
i=0; 
diff=toll+1; %Inizializzazione contatore

%% LOOP
while diff>=toll;
    
    %% Contatori
    toll=1e-100; %Errore concesso
    i=i+1;       %Contatore
    
    %% Calcolo variabili
    Tmedia=(Tin+Tout)/2; %Calcolo temperatura media
    cpi=calspec(Tmedia,c1,c2,c3,c4,c5,MWH2O); %Function per il calcolo del calore specifico
    Deltat=Q/(MH2O*cpi); %Calcolo differenza di temperatura 
    Tout=Tin+Deltat; %Calcolo temperatura in uscita
    
    %% Vettori
    x(i)=Tout; %Vettore iterate - Temperature in uscita
    z(i)=Deltat; %Vettore iterate - Differenza di temperatura
    y(i)=cpi; %Vettore iterate - Calore specifico
    
    %% Controllo
    if i>1
        diff=abs(x(i)-x(i-1)); %Calcolo dell'errore
    else
    end
        
      

end

%% DISPLAY
disp('TEMPERATURA IN USCITA');
disp(x(length(x)));
disp('ITERAZIONI A CONVERGENZA');
disp(length(x));
disp(' ');
%%disp('VETTORI');
%%disp(' ');
%%disp('Vettore temperature in uscita:');
%%disp(' ');
%%disp(x);
%%disp(' ');
%%disp('Vettore calori specifici:');
%%disp(' ');
%%disp(y);
%%disp(' ');
%%disp('Vettore differenze di temperatura:');
%%disp(' ');
%%disp(z);





                                    
                  


    
    
           
    
  
    
        
        
        
        
