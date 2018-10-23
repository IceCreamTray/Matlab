%% PROFILE OF PRESSURE IN A PUNCTURED TANK
% Finds the profile of pressure through transfer function and through analytical mass balance.
% ------------------------------------------------------------------------------------------

%% TRANSFER FUNCTION
% Constants
K = 2834000 % 1/ms
tau = 42.5 % s
w0 = 0.193 % Kg/s
P0 = (500 + 101.3)*10^3 % Pa

% Vectors
t = 0:300;

% Function
Pt = P0 - K .* (1 - exp(-t ./ tau)) .* w0

% Plot
plot(t,Pt*10^-3,'Linewidth',2);
title("Pressure profile");
ylabel("P(t) [kPa]");
xlabel("t [s]");

%% ACTUAL PROFILE
% Constants
V = 1.5; % m^3
A = 0.785 * 10^-4; % m^2
drho = 2 * 28.97 / 8.314 / 343.15 * 10^-3; % s^2/m^2
C = (drho / 2)^-1 / V;
Patm = 101300; % Pa

[Pf] = ode23(@odefun, [0 200], P0, A, drho, Patm, C);
disp(Pf)

function [dPdt] = odefun(t ,P, A, drho, Patm, C)
    dPdt = (- A * (drho * P * (P - Patm))^0.5 * C);
end
