% =========================== STEAM REFORMING =============================
%
% Calculates mass and energy balances for a packed bed real reactor.
% The molar flow rates are integrated upon the volume of catalyst.
% Plots molar flow rates profile along the reactor, molar fraction
% profile,temperature profile, heat flow of reaction.
% Provides graphical indications of methane conversion.
% It is possible to add a heat duty supply.
%
% The listed reactions are a simplified version of methane steam reforming:
% 1) CH4 + H2O <=> CO + 3H2 (methane steam reforming)
% 2) CO + H2O <=> CO2 + H2 (water-gas shift)
% 3) CH4 + 2H2O <=> CO2 + 4H2 (methanation)
% In this application, the overall density is not constant and the reactor 
% works at non isothermal conditions.
%
% =========================================================================


function Non_isothermal_reactor

	%% ===================== MASS AND ENERGY BALANCES =====================

	function [PFR, QR] = BMiBe(vol, Y, nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai)
		Ni = Y(1:5);							% kmol/h
		Tk = Y(6);

		% Define partial pressure as function of molar flowrates  
		Ntot=sum(Ni)
		p = P * Ni ./ Ntot						% i = CO,H2,CH4,H2O,CO2

		% Define kinetic constants
		kpj = [exp(30.481 - 27.187e+03 / Tk) exp(-3.924 + 4.291e+03 / Tk) exp(26.891 - 23.258e+03 / Tk)]';
		kj = Aj .* (exp(-Ej ./ (Rg * Tk)));
		ki = Ai .* (exp(-DHi ./ (Rg * Tk)));

		% Calculate reaction rates
		DEN = 1 + ki(1) * p(1) + ki(2) * p(2) + ki(3) * p(3) + ((ki(4) * p(4)) / p(2)); 
		R1 = ((kj(1) / (p(2) ^ 2.5)) * (p(3) * p(4) - ((p(2) ^ 3 * p(1)) / kpj(1)))) / ((DEN) ^ 2);
		R2 = ((kj(2) / p(2)) * (p(1) * p(4) - ((p(2) * p(5)) / kpj(2)))) / ((DEN) ^ 2);
		R3 = ((kj(3) / (p(2) ^ 3.5)) * ((p(3) * (p(4) ^ 2)) - (p(2) ^ 4 * p(5) / kpj(3)))) / ((DEN) ^ 2);
		R = [R1 R2 R3]';

		% Define entalphies of reaction
		DH1 = ((-4.47*(10^13)) * (Tk^ -4.459)) + 226.9;       
		DH2 = -271.4 * (Tk ^ -0.2977);
		DH3 = 99.52 * (Tk^0.0937);
		DHj = ([DH1 DH2 DH3])'					% Kj /mol    

		% Define specific heat
		cpi = ((ai.*Tk + bi).* Rg*10^-1)';
		cp_mix = sum(((Ni./Ntot)).*cpi);

		% Add heat duty if wanted
		Qdot = 300 * 36/19;

		% Adiabatic case
		%Qdot = 0;
		
		% Define heat of reaction
		QR = -R' * DHj;
		
		% Define total heat
		Q = Qdot + QR
		
		% Energy balance 
		Tfun = (Q * ro_bulk / phi) / cp_mix / Ntot;
		
		%Tfun = 0;								% For isothermal case

		% Define functions to integrate
		r = nu * R;
		PFR = [(r * (ro_bulk / phi)); Tfun];
		QR = abs(-R' * DHj);
		
	end

	%% =========================== MAIN SCRIPT ============================

	clc, clear all, close all;
	format long;

	% Initialize operating conditions
	T = 650 + 273.15;								% Kelvin
	P = 5;											% Bar

	% Initialize constants
	Rg = 8.3145 * 10^-2;							% Kj / K mol                                 

	% Initialize kinetic and equilibrium parameters
	Aj = [4.225e+15 1.955e+06 1.020e+15]';			% With j = 1,2,3 - bar^-1
	Ai = [8.23e-05 6.12e-09 6.65e-04 1.77e+05]';	% With i = CO,H2,CH4,H2O - bar^-1, adim
	Ej = [240.1 67.13 243.9]' * 10;					% Bar m^3/kmol
	DHi = [-70.65 -82.90 -38.28 88.68]' * 10;		% Bar m^3/kmol
	ro_bulk = 1900;									% Kg/m3_reactor
	phi = 0.6;										% Adim

	% Constants definition
	ai = [0.716 0.370 5.529 1.558 1.593] / 10^3;
	bi = [3.269 3.266 3.292 3.409 4.925];

	% Inlet molar flow rates
	Y0 = ([0 2.63.*3.6 50.*3.6 150.*3.6 0 T]);		% kmol/h, K         
	Ni0 = Y0(1:5);

	% Define stoichiometric matrix
	nu = [ 1 -1  0
		   3  1  4
		  -1  0 -1
		  -1 -1 -2
		   0  1  1];


	%% Integration
	[vol, Y] = ode15s(@BMiBe, [0 0.012], Y0, [], nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai);
	Ni = Y(:,1:5);
	Nitot = sum(Ni');
	Tk = Y(:,6);


	%% Plots

	%%
	% Define profile along the volume of catalyst
	Nprof = [Ni(:,1) Ni(:,2) Ni(:,3) Ni(:,4) Ni(:,5)];

	figure
	plot(vol * 10^3, Nprof, 'LineWidth', 1.5), hold on
	title('Flow rate profiles along the volume of catalyst');
	legend('CO', 'H2','H2O', 'CH4', 'CO2');
	xlabel('Volume of catalyst [l]');
	ylabel('Molar flow rate of species [kmol/h]');


	%% 
	% Total flowrate

	figure
	plot(vol * 10^3, Nitot, 'LineWidth', 1.5);
	title('Total flow rate profile along the volume of catalyst');
	xlabel('Volume of catalyst [l]');
	ylabel('Total molar flow rate [kmol/h]');


	%% 
	% Define molar fractions plot
	x = Ni ./ Nitot';

	figure 
	plot(vol * 10^3, x, 'Linewidth', 1.5), hold on
	xlabel('Volume of catalyst [l]');
	ylabel('Molar fractions of species');
	title('Molar fractions profiles along the volume of catalyst');
	legend('CO', 'H2','H2O', 'CH4', 'CO2');


	%%
	% Temprerature plot

	figure 
	plot(vol * 10^3, Tk - 273.15, 'LineWidth', 1.5)
	xlabel('Volume of catalyst [l]');
	ylabel('Temperature(°C)');
	title('Temperature profile along the volume of catalyst');


	%%
	% Define conversion plot
	conv_CH4 = (Ni0(3) - Ni(:,3)) ./ (Ni0(3));

	figure
	plot(vol * 10^3, conv_CH4, 'Linewidth', 1.5);
	l = refline(0, 0.5);
	set(l, 'LineStyle', '--', 'Color', 'black');
	title('Conversion of methane along reactor volume');
	xlabel('Volume of catalyst[l]');
	ylabel('Conversion of Methane');


	%% 
	% Heat of reaction
	for u = 1:length(vol)
		[PFR, QR] = BMiBe(vol(u), Y(u,:)', nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai);
		Qvec(u) = QR;
	end

	figure
	plot(vol * 10^3, Qvec, 'LineWidth', 1.5), hold on
	title('Heat of reaction vs volume of catalyst');
	xlabel('Volume of catalyst [l]');
	ylabel('Heat of reaction modulus [KJ]');
	xlim([0 2]);
	ylim([0 30000]);
	
	%%
	% Enthalpies temperature dependancy plot
	Ti = linspace(273,5000);
	for k = 1:length(Ti)
	DH1(k) = ((-4.47 * (10^13)) * (Ti(k)^ -4.459)) + 226.9;       
	DH2(k) = -271.4 * (Ti(k) ^ -0.2977);
	DH3(k) = 99.52 * (Ti(k)^0.0937);
	end
	
	figure
	plot(Ti - 273.15, DH1, 'LineWidth', 1.5), hold on
	plot(Ti - 273.15, DH2, 'LineWidth', 1.5), hold on
	plot(Ti - 273.15, DH3, 'LineWidth', 1.5), hold on;
	title('Dependancy of enthalpy upon temperature');
	xlabel('Temperature (°C)');
	ylabel('Enthalpy [KJ/mol]');
	xlim([0 4000]);
	
	%%
	% Boyfriend contribution
	lotsa = -10 : 0.01 : 10;
	lo = 16 * (power(sin(lotsa), 3));
	ve = (13 * cos(lotsa)) - (5 * cos(2 * lotsa)) - (2 * cos(3 * lotsa)) - (cos(4 * lotsa));
	
	figure
	plot(lo,ve);
	title('929');
	
end