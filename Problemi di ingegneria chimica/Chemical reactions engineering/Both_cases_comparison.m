%% ==================== MAIN SCRIPT ====================
format long;
clc, clear all, close all

% Initialize operating conditions
T = 650+273.15;                                   % Kelvin
P = 5;                                            % Pa

% Initialize constants
Rg = 8.3145*10^-3;                                  % KJ/ K mol

% Initialize kinetic and equilibrium parameters
Aj = [4.225e+15 1.955e+06 1.020e+15]';            % with j = 1,2,3 - bar^-1
Ai = [8.23e-05 6.12e-09 6.65e-04 1.77e+05]';      % with i = CO,H2,CH4,H2O - bar^-1, adim
Ej = [240.1 67.13 243.9]';                         % Kj /mol
DHi = [-70.65 -82.90 -38.28 88.68]';              % Kj/mol
ro_bulk = 1900;                                   % Kg/m3_reactor
phi = 0.6;                                        % adim

% Constants definition
ai = [0.716 0.370 5.529 1.558 1.593] / 10^3;
bi = [3.269 3.266 3.292 3.409 4.925];



%cpi calculations


% Calculate equilibrium constants, with j=1,2,3 - bar2, adim, bar2

% Calculate kinetic constants

% Inlet molar flow rates
Y0 = ([0*3.6 2.63.*3.6 50.*3.6 150.*3.6 0*3.6 T]);                % from mol/s to kmol/h           
Ni0 = Y0(1:5);
Niso0=Ni0;
% Define stoichiometric matrix
nu=[+1 -1 0
    +3 +1 +4
    -1 0 -1
    -1 -1 -2
    0 +1 +1];

%% Integration
%Non isothermal
[vol,Y] = ode15s(@BMiBe,[0 0.008],Y0,[],nu,phi,ro_bulk,P,Rg,ai,bi,DHi,Ej,Aj,Ai);
Ni = Y(:,1:5);
Nitot = sum(Ni');
Tk = Y(:,6);
%isothermal
[volu,Niso] = ode23(@BMiso,[0 0.008],Niso0,[],nu,phi,ro_bulk,P,Aj,Ai,Ej,DHi,T,Rg);

for u=1:length(vol)
	[PFR,QR] = BMiBe(vol(u),Y(u,1:6),nu,phi,ro_bulk,P,Rg,ai,bi,DHi,Ej,Aj,Ai);
	Qvec(u)=QR;
end


%% Plots
% Red is non iso, blue is iso.

%% Define profile along the volume of catalyst (B)
%non isothermal
figure

Nprof = [Ni(:,1) Ni(:,2) Ni(:,3) Ni(:,4) Ni(:,5)];
Nprofiso = [Niso(:,1) Niso(:,2) Niso(:,3) Niso(:,4) Niso(:,5)];
plot(vol*10^3, Nprof,'-r', 'LineWidth', 1.5), hold on
plot(volu*10^3, Nprofiso,'Color','Blue', 'LineWidth', 1.5);
xlabel('Volume of catalyst [L]');
ylabel('Molar flow rate of species [kmol/h]');
title('Flow rate profile along the volume of catalyst');



%% Total flowrate non iso
figure
plot(vol*10^3,Nitot,'-o','LineWidth',1.5);
xlabel('Volume of catalyst [L]');
ylabel('Total molar flow rate [kmol/h]');
title('Flow rate profile along the volume of catalyst');

%% Define molar fractions plot (A)
x = bsxfun(@rdivide,[Ni(:,1) Ni(:,2) Ni(:,3) Ni(:,4) Ni(:,5)],Nitot');
figure 
plot(vol*10^3,x,'-o','Linewidth',1.5), hold on
xlabel('Volume of catalyst [L]');
ylabel('Molar fractions of species');
title('Molar fractions profiles along the volume of catalyst');
legend('CO', 'H2','H2O', 'CH4', 'CO2');

%% Temprerature plot (C)
figure 
plot(vol*10^3,Tk-273.15,'-o','LineWidth',1.5)
xlabel('Volume of catalyst [L]');
ylabel('Temperature(°C)');
title('Temperature profile along the volume of catalyst');

%% Define conversion plot
figure
conv_CH4 = (Ni0(3) - Ni(:,3)) ./ (Ni0(3));
conviso_CH4 = (Niso0(3) - Niso(:,3)) ./ (Niso0(3));
plot(vol*10^3,conv_CH4,'-r','Linewidth',1.5),hold on
plot(volu*10^3, conviso_CH4,'Color','Blue', 'LineWidth', 1.5), hold on
xlabel('Volume of catalyst[L]');
ylabel('Conversion of Methane');
title('Conversion of methane along reactor volume');
l = refline(0, 0.5);
set(l, 'LineStyle', '--', 'Color', 'black');

figure
plot(Tk-273.15, conv_CH4, 'LineWidth', 1.5), hold on
xlabel('Temperature [C]');
ylabel('Conversion of Methane');
title('Conversion vs Temperature');









%% Heat of reaction
figure
plot(vol*10^3,Qvec*10^2,'-o','LineWidth',1.5);
xlabel('Volume of catalyst[L]');
ylabel('Heat of reaction [KJ]');
title('Heat of reaction vs volume of catalyst');
