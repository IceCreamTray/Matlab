function [PFR] = BMiso(volu,Niso,nu,phi,ro_bulk,P,Aj,Ai,Ej,DHi,T,Rg)
% Calculate equilibrium constants, with j=1,2,3 - bar2, adim, bar2
kpj = [exp(30.481 - 27.187e+03 / T) exp(-3.924 + 4.291e+03 / T) exp(26.891 - 23.258e+03 / T)]';

% Calculate kinetic constants
kj = Aj .* (exp(-Ej ./ (Rg * T)));
ki = Ai .* (exp(-DHi ./ (Rg * T)));

% Define partial pressure as function of molar flowrates  
Ntot=sum(Niso);
p = P * Niso ./ Ntot;       % i = CO,H2,CH4,H2O,CO2

% Calculate reaction rates
DEN = 1 + ki(1) * p(1) + ki(2) * p(2) + ki(3) * p(3) + ((ki(4) * p(4)) / p(2)); 
R1 = ((kj(1) / (p(2) ^ 2.5)) * (p(3) * p(4) - ((p(2) ^ 3 * p(1)) / kpj(1)))) / ((DEN) ^ 2);
R2 = ((kj(2) / p(2)) * (p(1) * p(4) - ((p(2) * p(5)) / kpj(2)))) / ((DEN) ^ 2);
R3 = ((kj(3) / (p(2) ^ 3.5)) * ((p(3) * (p(4) ^ 2)) - (p(2) ^ 4 * p(5) / kpj(3)))) / ((DEN) ^ 2);
R = [R1 R2 R3]';

% Calculate production rates
r = nu * R;

% Define function to integrate
PFR = r * (ro_bulk/phi);
end