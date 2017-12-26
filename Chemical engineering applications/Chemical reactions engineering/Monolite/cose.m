function  [time,Y]=cose(D,L,Sh,d,a,T,Creag,R,par,time)
% come CatNP_V0, ma per il sottocaso di F che si muove come PFR
% 
% hm è calcolato da correlazioni e sono parzialmente implementate diverse geometrie
% Cambiando correlazioni per Sh e formulazioni per a, può essere usato per
% reattori a letto impaccato, o a canale (es. monolita), o a rete
 
        %calcolo coeff di mass transfer 
  
		
       
        hm  =  Sh*D/d;  %m/s
    
   


% NB: la stima precedente può essere critica per D piccolo e/o reaz. veloce
CIN = Creag;

% risoluz. dei bilanci in fase fluida e all'interfaccia, lungo z, fino a L

M = [1 0
     0 0];          % ADE (ODE la 1 eq, AE la 2a)
options = odeset('Mass',M);
[time,Y] = ode15s(@BM,time,CIN,options,a,hm,T,par);

% disegna i profili adimensionali lungo z di Cbulk e Csuperficie (Cs)



function ADE = BM(time,y,a,hm,T,par) % =========================================
% formulazione dei bilanci in fase fluida e all'interfaccia
disp(y)
disp(a)
disp(hm)
disp(T)
disp(par)
Cb = y(1);
Cs = y(2);


% bilancio materiale in fase fluida (ODE)
% bilancio materiale in fase solida (all'interfaccia) (AE)

ADE = [    -a*hm * (Cb - Cs)
              hm * (Cb - Cs)  +   r(par,T,Cs)  ]
end
function r = r(par,T,Cs)  % ===================================================
% modello di velocità di reazione in funzione delle C in fase fluida,
% antistanti la superficie catalitica (Cs)
nu = -1;
% Vel. di reazione 
R    = k(par,T)*Cs;
r    = nu*R;
end

function k = k(par,T)   % ===================================================
% calcola la costante cinetica alla temperatura assegnata


A    = par(1);   %[mol/m^2*s]
Eatt = par(2);     %[J/mol]
R    = 8.314;   %[J/mol*K]
k    = A*exp(-Eatt/(R*T));
end
end

