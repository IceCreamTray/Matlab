%% ==================== STEAM REFORMING ====================
% Calculates mass balances for a packed bed differential reactor.
% The molar flow rates are integrated upon the volume of catalyst.
% Plots molar flow rates profile along the reactor,
% dry products conversion, conversion of methane, selectivities,
% yields, gas velocities, residence time.
% Provides graphical indications of catalyst volume, yields, selectivities, gas
% velocities, residence time at 50% methane conversion.

% The reactions occurring are a simplified version of steam reforming of methane:
% 1) CH4 + H2O <=> CO + 3H2 (methane steam reforming)
% 2) CO + H2O <=> CO2 + H2 (water-gas shift)
% 3) CH4 + 2H2O <=> CO2 + 4H2 (methanation)
% In this application, the overall density is not constant and the reactor 
% is considered to be working at isothermal conditions.

%% ==================== MAIN SCRIPT ====================
format long;
clc, clear all, close all

% Initialize operating conditions
T = 650+273.15;                                   % Kelvin
P = 5;                                            % bar

% Initialize constants
Rg = 8.3145e-03;                                  % bar/m3 K kmol

% Initialize kinetic and equilibrium parameters
Aj = [4.225e+15 1.955e+06 1.020e+15]';            % with j = 1,2,3 - bar^-1
Ai = [8.23e-05 6.12e-09 6.65e-04 1.77e+05]';      % with i = CO,H2,CH4,H2O - bar^-1, adim
Ej = [240.1 67.13 243.9]';                        % bar m3/kmol
DHi = [-70.65 -82.90 -38.28 88.68]';              % bar m3/kmol
ro_bulk = 1900;                                   % Kg/m3_reactor
phi = 0.6;                                        % adim

% Calculate equilibrium constants, with j=1,2,3 - bar2, adim, bar2
kpj = [exp(30.481 - 27.187e+03 / T) exp(-3.924 + 4.291e+03 / T) exp(26.891 - 23.258e+03 / T)]';

% Calculate kinetic constants
kj = Aj .* (exp(-Ej ./ (Rg * T)));
ki = Ai .* (exp(-DHi ./ (Rg * T)));

% Inlet molar flow rates
Ni0 = ([0 2.63 50 150 0]' .* 3.6);                % from mol/s to kmol/h           

% Define stoichiometric matrix
nu=[+1 -1 0
    +3 +1 +4
    -1 0 -1
    -1 -1 -2
    0 +1 +1];

%% Integration
[vol,Ni] = ode23(@BMi,[0 0.008],Ni0,[],ki,kj,kpj,nu,phi,ro_bulk,P);

%% Plots

% Define profile along the volume of catalyst
figure(1)

Nprof = [Ni(:,1) Ni(:,2) Ni(:,3) Ni(:,4) Ni(:,5)];
plot(vol, Nprof, '-o', 'LineWidth', 1.5), hold on
legend('CO', 'H2','H2O', 'CH4', 'CO2');
xlabel('Volume of catalyst [m^3]');
ylabel('Molar flow rate of species [kmol/h]');
title('Flow rate profile along the volume of catalyst');

%%
% Define dry products plot
Ndry = [Ni(:,1), Ni(:,2), Ni(:,3), Ni(:,5)]';
Ndry_tot = sum(Ndry);
xdry = bsxfun(@rdivide, Ndry, Ndry_tot);

figure(2)

plot(vol, xdry, '-o', 'LineWidth', 1.5), hold on
xlabel('Volume of catalyst [m^3]');
ylabel('Fraction of dry products');
title('Dry products profile along the volume of catalyst');
legend('CO', 'H2', 'CH4', 'CO2');

%%
% Define conversion plot
conv_CH4 = (Ni0(3) - Ni(:,3)) ./ (Ni0(3));

figure(3)

plot(vol, conv_CH4, 'LineWidth', 1.5), hold on
xlabel('Volume of catalyst [m^3]');
ylabel('Conversion of Methane');
title('Conversion of Methane along the volume of catalyst');
l = refline(0, 0.5);
set(l, 'LineStyle', '--', 'Color', 'black');

%%
% Define yields plots
yi_h2_ch4 = ((Ni(:,2) - Ni0(2)) ./ (Ni0(3))) .* (1/3.5);
yi_co_ch4 = (Ni(:,1) ./ Ni0(3));

figure (4)

plot (vol, yi_co_ch4, 'LineWidth', 1.5), hold on, grid on
plot(vol, yi_h2_ch4, '-r', 'LineWidth', 1.5)
l = line([4.18e-03 4.18e-03], [0 1]);
set(l, 'LineStyle', '--', 'Color', 'black');
xlabel('Volume of catalyst [m^3]');
ylabel('Yield of H2 and CO on CH4');
title('Yield of H2 and CO on CH4 along the volume of catalyst');
legend('CO','H2');


%%
%Define selectivities plot
sh2 = (Ni(:,2) - Ni0(2)) / (Ni0(3) - Ni(:,3)) * 1/4;
sco = (Ni(:,1) - Ni0(1)) / (Ni0(3) - Ni(:,3));

figure (5)

plot(vol, sh2(:,46), '-r', 'LineWidth', 1.5), hold on, grid on
plot(vol, sco(:,46), 'LineWidth', 1.5) 
legend('H2','CO');
l = line([4.18e-03 4.18e-03], [0 1]);
set(l, 'LineStyle', '--', 'Color', 'black');
xlabel('Volume of catalyst [m^3]');
ylabel('Selectivity of H2,CO on CH4');
title('Selectivities of H2,CO on CH4 along the volume of catalyst');


%%
% Ratio between inlet and outlet gas velocity
vo = (Ni.* (Rg * T)) ./ P;
volflow_out_tot = vo(:,1) + vo(:,2) + vo(:,3) + vo(:,4) + vo(:,5);
vi = (Ni0 .* (Rg * T)) ./ P;
volflow_in_tot = vi(1) + vi(2) + vi(3) + vi(4) + vi(5);
ratio = (volflow_out_tot / volflow_in_tot);

figure(6)

plot(vol, ratio, 'LineWidth', 1.5), hold on, grid on
l = line([4.18e-03 4.18e-03], [1 1.5]);
set(l, 'LineStyle', '--', 'Color', 'black');
xlabel('Volume of catalyst [m^3]');
ylabel('Outlet gas velocity and inlet gas velocity ratio');
title('Outlet and inlet gas velocity ratio along the volume of catalyst');

%%
% Residence time
volu = ((Ni * 3.6).* (Rg ./ 1000 * T)) ./ P;
volflow_tot = volu(:,1) + volu(:,2) + volu(:,3) + volu(:,4) + volu(:,5);
h = cumtrapz(vol,1 ./ (volflow_tot));

figure(7)

plot(vol, h * 10 ^ 3,'-o', 'LineWidth', 1.5), grid on
l = line([4.18e-03 4.18e-03],[0 2]);
set(l, 'LineStyle', '--', 'Color', 'black');
xlabel('Volume of catalyst [m^3]');
ylabel('Residence time [ms]');
title('Residence time profile along the volume of catalyst');

tau = trapz(vol, 1 ./ volflow_tot);
format short
disp('Residence time value at output [ms]: ');
disp(tau * 10 ^ 3);
