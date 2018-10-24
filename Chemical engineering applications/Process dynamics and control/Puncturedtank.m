%% PROFILE OF PRESSURE IN A PUNCTURED TANK
% Finds the profile of pressure through transfer function and through analytical mass balance.
% ------------------------------------------------------------------------------------------

clc, clear all, close all;

%% PROFILE APPROXIMATION THROUGH TRANSFER FUNCTION
% Constants
K = 2834000; % 1/ms
tau = 42.5; % s
w0 = 0.193; % Kg/s
P0 = (500 + 101.3)*10^3; % Pa

% Vectors
t1 = 0:300;

% Function in time domain
Pt = P0 - K .* (1 - exp(-t1 ./ tau)) .* w0;


%% ACTUAL PROFILE
% Constants
V = 1.5; % m^3
A = 0.785 * 10^-4; % m^2
drho = 2 * 28.97 / 8.314 / 343.15 * 10^-3; % s^2/m^2
C = (drho / 2)^-1 / V;
Patm = 101300; % Pa

% Material balance solution
[t,P] = ode23(@Bmi, [0 300], P0, [], A, drho, Patm,V);

% TIME THE TANK TAKES TO GET EMPTY
% Calculated through the integral of the material balance
f = @(P) (drho ./ 2 .* V)./(-A.*(drho.*P.*(P-Patm)).^0.5)
res = integral(f,601300,101300)
% Calculated through the reverse transform
tr = (log(1-((P0 - Patm)/w0/K))*tau)*-1


%% PLOTS
plot(t1,Pt*10^-3,'Linewidth',2),  hold on
title("Pressure profile");
ylabel("P [kPa]");
xlabel("t [s]");
plot(t,P*10^-3,'Linewidth', 2);
l = refline(0, 101.3);
set(l, 'LineStyle', '--', 'Color', 'black');
legend('Approximated profile','Actual profile', 'Atmospheric pressure line');


function [dPdt] = Bmi (t,P,A,drho,Patm,V)
    dPdt = (-A*(drho*P*(P-Patm))^0.5)*(drho/2)^-1*V^-1;
end
