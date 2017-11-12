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

% Constants definition
ai = [0.716 0.370 5.529 1.558 1.593] * 10^3;
bi = [3.269 3.266 3.292 3.409 4.925];



%cpi calculations


% Calculate equilibrium constants, with j=1,2,3 - bar2, adim, bar2

% Calculate kinetic constants

% Inlet molar flow rates
Y0 = ([0 2.63.*3.6 50.*3.6 150.*3.6 0 T]');                % from mol/s to kmol/h           
Ni0 = Y0(1:5);
% Define stoichiometric matrix
nu=[+1 -1 0
    +3 +1 +4
    -1 0 -1
    -1 -1 -2
    0 +1 +1];

%% Integration
[vol,Y] = ode15s(@BMiBe,[0 0.008],Y0,[],nu,phi,ro_bulk,P,Rg,ai,bi,DHi,Ej,Aj,Ai);
Ni = Y(:,1:5);
Nitot = sum(Ni');
Tk = Y(:,6);


%% Plots

% Define profile along the volume of catalyst (B)
figure

Nprof = [Ni(:,1) Ni(:,2) Ni(:,3) Ni(:,4) Ni(:,5)];
plot(vol, Nprof, '-o', 'LineWidth', 1.5), hold on
legend('CO', 'H2','H2O', 'CH4', 'CO2');
xlabel('Volume of catalyst [m^3]');
ylabel('Molar flow rate of species [kmol/h]');
title('Flow rate profile along the volume of catalyst');

figure
plot(vol,Nitot,'-o','LineWidth',1.5);
xlabel('Volume of catalyst [m^3]');
ylabel('Total molar flow rate [kmol/h]');
title('Flow rate profile along the volume of catalyst');

% Define molar fractions plot (A)
x = bsxfun(@rdivide,[Ni(:,1) Ni(:,2) Ni(:,3) Ni(:,4) Ni(:,5)],Nitot');
figure 
plot(vol,x,'-o','Linewidth',1.5), hold on
xlabel('Volume of catalyst [m^3]');
ylabel('Molar fractions of species');
title('Molar fractions profiles along the volume of catalyst');
legend('CO', 'H2','H2O', 'CH4', 'CO2');

% Temprerature plot (C)
figure 
plot(vol,Tk-273.15,'-o','LineWidth',1.5)
xlabel('Volume of catalyst [m^3]');
ylabel('Temperature(°C)');
title('Temperature profile along the volume of catalyst');
adiab_delta = (Tk - (650+273.15));
figure
plot(vol,adiab_delta,'-o','LineWidth',1.5);
xlabel('Volume of catalyst [m^3]');
ylabel('Temperature change');
title('Temperature change along the volume of catalyst');


% %%
% % Define dry products plot
% Ndry = [Ni(:,1), Ni(:,2), Ni(:,3), Ni(:,5)]';
% Ndry_tot = sum(Ndry);
% xdry = bsxfun(@rdivide, Ndry, Ndry_tot);
% 
% figure(2)
% 
% plot(vol, xdry, '-o', 'LineWidth', 1.5), hold on
% xlabel('Volume of catalyst [m^3]');
% ylabel('Fraction of dry products');
% title('Dry products profile along the volume of catalyst');
% legend('CO', 'H2', 'CH4', 'CO2');
% 
% %%
% Define conversion plot
figure
conv_CH4 = (Ni0(3) - Ni(:,3)) ./ (Ni0(3));
plot(Tk-273.15, conv_CH4, 'LineWidth', 1.5), hold on
xlabel('Temperature [C]');
ylabel('Conversion of Methane');
title('Conversion vs Temperature');

figure
plot(vol,conv_CH4,'Linewidth',1.5);
xlabel('Volume of catalyst[m^3]');
ylabel('Conversion of Methane');
title('Conversion of methane along reactor volume');
l = refline(0, 0.5);
set(l, 'LineStyle', '--', 'Color', 'black');



% 
% %%
% % Define yields plots
% yi_h2_ch4 = ((Ni(:,2) - Ni0(2)) ./ (Ni0(3))) .* (1/3.5);
% yi_co_ch4 = (Ni(:,1) ./ Ni0(3));
% 
% figure (4)
% 
% plot (vol, yi_co_ch4, 'LineWidth', 1.5), hold on, grid on
% plot(vol, yi_h2_ch4, '-r', 'LineWidth', 1.5)
% l = line([4.18e-03 4.18e-03], [0 1]);
% set(l, 'LineStyle', '--', 'Color', 'black');
% xlabel('Volume of catalyst [m^3]');
% ylabel('Yield of H2 and CO on CH4');
% title('Yield of H2 and CO on CH4 along the volume of catalyst');
% legend('CO','H2');
% 
% 
% %%
% %Define selectivities plot
% sh2 = (Ni(:,2) - Ni0(2)) ./ (Ni0(3) - Ni(:,3)) * 1/4;
% sco = (Ni(:,1) - Ni0(1)) ./ (Ni0(3) - Ni(:,3));
% 
% figure (5)
% 
% plot(vol, sh2, '-r', 'LineWidth', 1.5), hold on, grid on
% plot(vol, sco, 'LineWidth', 1.5) 
% legend('H2','CO');
% l = line([4.18e-03 4.18e-03], [0 1]);
% set(l, 'LineStyle', '--', 'Color', 'black');
% xlabel('Volume of catalyst [m^3]');
% ylabel('Selectivity of H2,CO on CH4');
% title('Selectivities of H2,CO on CH4 along the volume of catalyst');
% 
% 
% %%
% % Ratio between inlet and outlet gas velocity
% vo = (Ni.* (Rg * T)) ./ P;
% volflow_out_tot = sum(vo,2);
% vi = (Ni0 .* (Rg * T)) ./ P;
% volflow_in_tot = vi(1) + vi(2) + vi(3) + vi(4) + vi(5);
% ratio = (volflow_out_tot / volflow_in_tot);
% 
% figure(6)
% 
% plot(vol, ratio, 'LineWidth', 1.5), hold on, grid on
% l = line([4.18e-03 4.18e-03], [1 1.5]);
% set(l, 'LineStyle', '--', 'Color', 'black');
% xlabel('Volume of catalyst [m^3]');
% ylabel('Outlet gas velocity and inlet gas velocity ratio');
% title('Outlet and inlet gas velocity ratio along the volume of catalyst');
% 
% %%
% % Residence time
% volu = ((Ni.* 3.6).* ((Rg/1000) * T)) ./ P;
% volflow_tot = sum(volu,2);
% h = cumtrapz(vol,1 ./ (volflow_tot));
% 
% figure(7)
% 
% plot(vol, h * 10^3,'-o', 'LineWidth', 1.5), grid on
% l = line([4.18e-03 4.18e-03],[0 2]);
% set(l, 'LineStyle', '--', 'Color', 'black');
% xlabel('Volume of catalyst [m^3]');
% ylabel('Residence time [ms]');
% title('Residence time profile along the volume of catalyst');
% 
% for count=1:length(vol)
%     if vol(count) < 4.20e-03
%     vol2(count) = vol(count);
%     volflow2(count) = volflow_tot(count);
%     end
% end
% 
% tau = trapz(vol2, 1 ./ volflow2);
% format short
% disp('Residence time value at output [ms]: ');
% disp(tau * 10 ^ 3);
