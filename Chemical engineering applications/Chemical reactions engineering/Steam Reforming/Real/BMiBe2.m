function [PFR,QR] = BMiBe(vol,Y,nu,phi,ro_bulk,P,Rg,ai,bi,DHi,Ej,Aj,Ai)
Ni = Y(1:5);										% kmol/h
disp('Ni'),disp(Ni);
Tk = Y(6);
%disp('Tk'), disp(Tk)% K
% Define partial pressure as function of molar flowrates  
Ntot=sum(Ni)
p = P * Ni ./ Ntot       % i = CO,H2,CH4,H2O,CO2


kpj = [exp(30.481 - 27.187e+03 / Tk) exp(-3.924 + 4.291e+03 / Tk) exp(26.891 - 23.258e+03 / Tk)]';
kj = Aj .* (exp(-Ej ./ (Rg * Tk)));
ki = Ai .* (exp(-DHi ./ (Rg * Tk)));

% Calculate reaction rates
DEN = 1 + ki(1) * p(1) + ki(2) * p(2) + ki(3) * p(3) + ((ki(4) * p(4)) / p(2)); 
R1 = ((kj(1) / (p(2) ^ 2.5)) * (p(3) * p(4) - ((p(2) ^ 3 * p(1)) / kpj(1)))) / ((DEN) ^ 2);
R2 = ((kj(2) / p(2)) * (p(1) * p(4) - ((p(2) * p(5)) / kpj(2)))) / ((DEN) ^ 2);
R3 = ((kj(3) / (p(2) ^ 3.5)) * ((p(3) * (p(4) ^ 2)) - (p(2) ^ 4 * p(5) / kpj(3)))) / ((DEN) ^ 2);
R = [R1 R2 R3]';


DH1 = ((-4.47*(10^13)) * (Tk^ -4.459)) + 226.9;       
DH2 = -271.4 * (Tk ^ -0.2977);
DH3 = 99.52 * (Tk^0.0937);
DHj = ([DH1 DH2 DH3])'.*10  							% bar m^3 /kmol            
cpi = ((ai.*Tk + bi).* Rg)'

             
vtmp = ((Ni./Ntot).*cpi)
cp_mix = sum(vtmp);


%Q = 300*10   %KW/L - Kj/s*dm^3 -   bar m^3/s*m^3                    
Q=0;

% Energy balance 
Tfun =0;

% Define function to integrate
r = nu * R
rlen = length(r);
PFR = [(r * (ro_bulk/phi)); Tfun]; %
QR = R'*DHj;                               
%QR = Q;
end