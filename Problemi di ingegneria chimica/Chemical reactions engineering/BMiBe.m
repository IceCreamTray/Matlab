function [PFR] = BMiBe(vol,Y,nu,phi,ro_bulk,P,Rg,ai,bi,DHi,Ej,Aj,Ai)
Ni = Y(1:5);
Tk = Y(6);
kpj = [exp(30.481 - 27.187e+03 / Tk) exp(-3.924 + 4.291e+03 / Tk) exp(26.891 - 23.258e+03 / Tk)]';
kj = Aj .* (exp(-Ej ./ (Rg .* Tk)));
ki = Ai .* (exp(-DHi ./ (Rg .* Tk)));

% Define partial pressure as function of molar flowrates  
Ntot=sum(Ni);
p = P * Ni ./ Ntot;       % i = CO,H2,CH4,H2O,CO2

% Calculate reaction rates
DEN = 1 + ki(1) * p(1) + ki(2) * p(2) + ki(3) * p(3) + ((ki(4) * p(4)) / p(2)); 
R1 = ((kj(1) / (p(2) ^ 2.5)) * (p(3) * p(4) - ((p(2) ^ 3 * p(1)) / kpj(1)))) / ((DEN) ^ 2);
R2 = ((kj(2) / p(2)) * (p(1) * p(4) - ((p(2) * p(5)) / kpj(2)))) / ((DEN) ^ 2);
R3 = ((kj(3) / (p(2) ^ 3.5)) * ((p(3) * (p(4) ^ 2)) - (p(2) ^ 4 * p(5) / kpj(3)))) / ((DEN) ^ 2);
R = [R1 R2 R3]';

% Calculate production rates
r = nu * R;

% define mole fracts
DH1 = ((-4.47*(10^13)) * (Tk^ -4.459)) + 226.9;
DH2 = -271.4 * (Tk^ -0.2977);
DH3 = 99.52 * (Tk^0.0937);
DHj = ([DH1 DH2 DH3]*10^5)' ;
cpi = (((ai.*Tk + bi).* Rg))';
xi = bsxfun(@rdivide,Ni,Ntot);
counter = bsxfun(@rdivide,xi,1./cpi);
cp_mix = sum(counter);
somma = bsxfun(@rdivide,R,1./DHj);
QR = sum(somma);
Q = 0;
%Energy balance 
contribute1 = ((Q*phi) + QR*vol);
contribute2 = (Ntot * Rg * Tk) / (P * cp_mix * phi);
Tfun =(contribute1*contribute2);

% Define function to integrate
PFR = [(r * (ro_bulk/phi)); Tfun];
end