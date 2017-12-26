function CatNP_PFR
% come CatNP_V0, ma per il sottocaso di F che si muove come PFR
% 
% hm è calcolato da correlazioni e sono parzialmente implementate diverse geometrie
% Cambiando correlazioni per Sh e formulazioni per a, può essere usato per
% reattori a letto impaccato, o a canale (es. monolita), o a rete
 
clc,clear all
close all

% dati da gas
D   = 1e-5;   % m^2/s
% v   = 1;      % m/s
% rho = 1;    % kg/m3
% mu  = 1e-5; % Pa*s
% 
% % dati da liq
% D = 1e-9;     % m^2/s
% v = 0.1;      % m/s
% rho = 1e3;    % kg/m3
% mu =  1e-3;   % Pa*s

%config = 1;     % letto impaccato
config = 2;     % canale con cat sulla sup. interna
% config = 3;     % rete (incompleta)

L = 8.5e-03;      % lunghezza reattore, m
T =700;    % K
%Sc  =  mu/(rho*D);

switch config
%     case 1           % letto impaccato
%         
%         dp = 2e-3;       % m, diametro pellet
%         dt = 2.5e-2;     % m, diametro tubo
%         eb = 0.35;       % grado di vuoto
%         a = 6/dp*(1-eb)/eb;  % 1/m  area specifica
% 
%         %calcolo coeff di mass transfer
%         vs = v;          % vel superificiale = v assunta sopra
%         Re  = (rho*vs*dp)/mu;
%     %     Jd  =  1/eb*(0.765/Re^0.82  + 0.365/Re^0.386);       % corr 11-69 di Fogler
%         Jd  =  2.19/Re^(2/3) + 0.78/Re^0.381;       % corr 14.5-2 Bird
%         Sh  =  Jd*Sc^(1/3)*Re;
%         hm  =  Sh*D*(1-eb)/dp;  %m/s

    case 2          % canale con cat sulla sup. interna
        
        d = 8e-3;       % m, diametro canale
        a = 4/d;        % 1/m  area specifica
        
        %calcolo coeff di mass transfer 
%         Re  = (rho*v*d)/mu;
%         if Re < 2100    % laminare
%             Sh  =  1.62  * (Re*Sc * d/L)^(1/3);
%         else
%             Sh  =  0.026 * Re^0.8 *Sc^(1/3);
%         end
Sh = 3;
        hm  =  Sh*D/d;  %m/s
    
    case 3          % rete ;         da costruire
end

CbIN = 4;        % conc di bulk in z=0, moli/m3
Cs   = CbIN;        % stima conc alla superficie in z=0, moli/m3
% NB: la stima precedente può essere critica per D piccolo e/o reaz. veloce
CIN =[CbIN Cs]';

% risoluz. dei bilanci in fase fluida e all'interfaccia, lungo z, fino a L
disp(sprintf('k"/hm = %g (possibly inconsistent units!)',k(T)/hm))
M = [1 0
     0 0];          % ADE (ODE la 1 eq, AE la 2a)
options = odeset('Mass',M);
[z,Y] = ode15s(@BM,[0 200],CIN,options,a,hm,T);
disp(Y)
% disegna i profili adimensionali lungo z di Cbulk e Csuperficie (Cs)
plot(z,Y)
xlabel('z [m]'), ylabel('C reagente/Cin'),legend('bulk','sup')



function ADE = BM(z,y,a,hm,T) % =========================================
% formulazione dei bilanci in fase fluida e all'interfaccia

Cb = y(1);
Cs = y(2);

% bilancio materiale in fase fluida (ODE)
% bilancio materiale in fase solida (all'interfaccia) (AE)

ADE = [    -a*hm * (Cb - Cs)
              hm * (Cb - Cs)   +   r(Cs,T)  ]

function r = r(Cs,T)  % ===================================================
% modello di velocità di reazione in funzione delle C in fase fluida,
% antistanti la superficie catalitica (Cs)

% Vel. di reazione 
n    = 1;
R    = k(T)*Cs^n;
r    = -n*R;

function k = k(T)   % ===================================================
% calcola la costante cinetica alla temperatura assegnata


A    = 1e+10;   %[mol/m^2*s]
Eatt = 1e5;     %[J/mol]
R    = 8.314;   %[J/mol*K]
k    = A*exp(-Eatt/(R*T));