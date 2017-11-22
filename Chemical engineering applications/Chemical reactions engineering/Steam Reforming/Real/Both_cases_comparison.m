% =========================== STEAM REFORMING =============================
%
% Calculates mass and energy balances for a packed bed real reactor.
% The molar flow rates are integrated upon the volume of catalyst.
% Plots molar flow rates profile along the reactor, conversions, yields,
% selectivities, heat of reaction, residence time.
% Provides graphical indications of methane conversion.
% Compares an adiabatic reactor (red) with an isotermal one (blue) and one 
% at constant heat supply (green).
% The listed reactions are a simplified version of methane steam reforming:
% 1) CH4 + H2O <=> CO + 3H2 (methane steam reforming)
% 2) CO + H2O <=> CO2 + H2 (water-gas shift)
% 3) CH4 + 2H2O <=> CO2 + 4H2 (methanation)
% In this application, the overall density is not constant.
%
% =========================================================================

function Comparison
	%% ========================= ADIABATIC CASE ===========================
	
	function[PFR, QR] = Adiabatic(vol, Y, nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai)
		Ni = Y(1:5);							% Kmol/h
		Tk = Y(6);

		% Define partial pressure as function of molar flowrates  
		Ntot = sum(Ni);
		p = P * Ni ./ Ntot;						% i = CO,H2,CH4,H2O,CO2

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
		DHj = [DH1 DH2 DH3]';					% Kj /mol    

		% Define specific heat
		cpi = ((ai .* Tk + bi).* Rg * 10^-1)';
		cp_mix = sum(((Ni ./ Ntot)) .* cpi);

		% Adiabatic case
		Qdot = 0;

		% Define heat of reaction
		QR = -R' * DHj;

		% Define total heat
		Q = Qdot + QR;

		% Energy balance 
		Tfun = (Q * ro_bulk / phi) / cp_mix / Ntot;

		% Define functions to integrate
		r = nu * R;
		PFR = [(r * (ro_bulk / phi)); Tfun];
		
		% Define heat of reaction
		QR = abs(-R' * DHj);

	end

	%% ========================= ISOTHERMAL CASE ==========================

	function[PFR2, QR2] = Isothermal(volu, Y, nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai)
		Niso = Y(1:5);							% Kmol/h
		Tkiso = Y(6);

		% Define partial pressure as function of molar flowrates  
		Ntotiso = sum(Niso);
		p = P * Niso ./ Ntotiso;				% i = CO,H2,CH4,H2O,CO2

		% Define kinetic constants
		kpj = [exp(30.481 - 27.187e+03 / Tkiso) exp(-3.924 + 4.291e+03 / Tkiso) exp(26.891 - 23.258e+03 / Tkiso)]';
		kj = Aj .* (exp(-Ej ./ (Rg * Tkiso)));
		ki = Ai .* (exp(-DHi ./ (Rg * Tkiso)));

		% Calculate reaction rates
		DEN = 1 + ki(1) * p(1) + ki(2) * p(2) + ki(3) * p(3) + ((ki(4) * p(4)) / p(2)); 
		R1 = ((kj(1) / (p(2) ^ 2.5)) * (p(3) * p(4) - ((p(2) ^ 3 * p(1)) / kpj(1)))) / ((DEN) ^ 2);
		R2 = ((kj(2) / p(2)) * (p(1) * p(4) - ((p(2) * p(5)) / kpj(2)))) / ((DEN) ^ 2);
		R3 = ((kj(3) / (p(2) ^ 3.5)) * ((p(3) * (p(4) ^ 2)) - (p(2) ^ 4 * p(5) / kpj(3)))) / ((DEN) ^ 2);
		R = [R1 R2 R3]';

		% Define entalphies of reaction
		DH1 = ((-4.47*(10^13)) * (Tkiso^ -4.459)) + 226.9;       
		DH2 = -271.4 * (Tkiso ^ -0.2977);
		DH3 = 99.52 * (Tkiso ^ 0.0937);
		DHj = [DH1 DH2 DH3]';					% Kj /mol    

		% Define specific heat
		cpi = ((ai .* Tkiso + bi).* Rg * 10^-1)';
		cp_mix = sum(((Niso ./ Ntotiso)) .* cpi);


		% Define heat of reaction
		QR2 = -R' * DHj;

		% Isothermal case
		Qdot = -QR2;
		
		% Define total heat
		Q = Qdot + QR2;

		% Energy balance 
		Tfun = (Q * ro_bulk / phi) / cp_mix / Ntotiso;

		% Define functions to integrate
		r = nu * R;
		PFR2 = [(r * (ro_bulk / phi)); Tfun];
		
		% Define heat of reaction
		QR2 = abs(-R' * DHj);

	end

	%% ========================= HEAT FLOW CONST ==========================
	
	function [PFR3, QR3] = HeatFlux(vola, Y, nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai)
		Nif = Y(1:5);								% Kmol/h
		Tkf = Y(6);

		% Define partial pressure as function of molar flowrates  
		Ntotf = sum(Nif);
		p = P * Nif ./ Ntotf;						% i = CO,H2,CH4,H2O,CO2

		% Define kinetic constants
		kpj = [exp(30.481 - 27.187e+03 / Tkf) exp(-3.924 + 4.291e+03 / Tkf) exp(26.891 - 23.258e+03 / Tkf)]';
		kj = Aj .* (exp(-Ej ./ (Rg * Tkf)));
		ki = Ai .* (exp(-DHi ./ (Rg * Tkf)));

		% Calculate reaction rates
		DEN = 1 + ki(1) * p(1) + ki(2) * p(2) + ki(3) * p(3) + ((ki(4) * p(4)) / p(2)); 
		R1 = ((kj(1) / (p(2) ^ 2.5)) * (p(3) * p(4) - ((p(2) ^ 3 * p(1)) / kpj(1)))) / ((DEN) ^ 2);
		R2 = ((kj(2) / p(2)) * (p(1) * p(4) - ((p(2) * p(5)) / kpj(2)))) / ((DEN) ^ 2);
		R3 = ((kj(3) / (p(2) ^ 3.5)) * ((p(3) * (p(4) ^ 2)) - (p(2) ^ 4 * p(5) / kpj(3)))) / ((DEN) ^ 2);
		R = [R1 R2 R3]';

		% Define entalphies of reaction
		DH1 = ((-4.47*(10^13)) * (Tkf^ -4.459)) + 226.9;       
		DH2 = -271.4 * (Tkf ^ -0.2977);
		DH3 = 99.52 * (Tkf ^ 0.0937);
		DHj = [DH1 DH2 DH3]';						% Kj /mol    

		% Define specific heat
		cpi = ((ai .* Tkf + bi).* Rg * 10^-1)';
		cp_mix = sum(((Nif ./ Ntotf)) .* cpi);

		% Add heat duty if wanted
		Qdot = 300 * 36/19;

		% Define heat of reaction
		QR3 = -R' * DHj;

		% Define total heat
		Q = Qdot + QR3;

		% Energy balance 
		Tfun = (Q * ro_bulk / phi) / cp_mix / Ntotf;


		% Define functions to integrate
		r = nu * R;
		PFR3 = [(r * (ro_bulk / phi)); Tfun];
		
		% Define heat of reaction
		QR3 = abs(-R' * DHj);

	end

	%% =========================== MAIN SCRIPT ============================
	
	clc, clear all, close all;
	format long;

	% Initialize operating conditions
	T = 650 + 273.15;									% Kelvin
	P = 5;												% Bar

	% Initialize constants
	Rg = 8.3145 * 10^-2;								% M^3 bar/ K kmol

	% Initialize kinetic and equilibrium parameters
	Aj = [4.225e+15 1.955e+06 1.020e+15]';				% With j = 1,2,3 - bar^-1
	Ai = [8.23e-05 6.12e-09 6.65e-04 1.77e+05]';		% With i = CO,H2,CH4,H2O - bar^-1, adim
	Ej = [240.1 67.13 243.9]' * 10;						% Bar m^3 /kmol
	DHi = [-70.65 -82.90 -38.28 88.68]' * 10;			% Bar m^3/mol
	ro_bulk = 1900;										% Kg/m3_reactor
	phi = 0.6;											% Adim

	% Constants definition
	ai = [0.716 0.370 5.529 1.558 1.593] / 10^3;
	bi = [3.269 3.266 3.292 3.409 4.925];

	% Inlet molar flow rates
	Y0 = ([0 2.63.*3.6 50.*3.6 150.*3.6 0 T]);			% Kmol/h           
	Ni0 = Y0(1:5);
	Niso0 = Ni0;

	% Define stoichiometric matrix
	nu = [+1 -1  0
		   3  1  4
		  -1  0 -1
		  -1 -1 -2
		   0  1  1];

	%% Integration

	% Non isothermal
	[vol, Y] = ode15s(@Adiabatic, [0 0.030], Y0, [], nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai);
	Ni = Y(:, 1:5);
	Nitot = sum(Ni');
	Tk = Y(:, 6);

	% Isothermal
	[volu, Y2] = ode15s(@Isothermal, [0 0.030], Y0, [], nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai);
	Niso = Y2(:, 1:5);
	Ntotiso = sum(Niso');
	Tkiso = Y2(:, 6);

	% Heat Flux
	[vola, Y3] = ode15s(@HeatFlux, [0 0.030], Y0, [], nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai);
	Nif = Y3(:, 1:5);
	Ntotf = sum(Nif');
	Tkf = Y3(:, 6);

	%% Plots
	% Red for adiabatic, blue for isothermal, green for constant heat supply.

	%%
	% Profile along the volume of catalyst
	Nprof = [Ni(:, 1) Ni(:, 2) Ni(:, 3) Ni(:, 4) Ni(:, 5)];
	Nprofiso = [Niso(:, 1) Niso(:, 2) Niso(:, 3) Niso(:, 4) Niso(:, 5)];
	%Nproff = [Nif(:,1) Nif(:,2) Nif(:,3) Nif(:,4) Nif(:,5)];

	figure
	plot(vol * 10^3, Nprof, '-r', 'LineWidth', 1.5), hold on
	plot(volu * 10^3, Nprofiso, 'Color', 'blue', 'LineWidth', 1.5), hold on
	%plot(vola * 10^3, Nproff, '-g', 'LineWidth', 1.5), hold on
	title('Flow rate profile along the volume of catalyst');
	xlabel('Volume of catalyst [l]');
	ylabel('Molar flow rate of species [kmol/h]');

	%% 
	% Define dry molar fractions plot
	x = bsxfun(@rdivide, [Ni(:, 1) Ni(:, 2) Ni(:, 3) Ni(:, 5)], Nitot');
	Nisodry = [Niso(:, 1), Niso(:, 2), Niso(:, 3), Niso(:, 5)]';
	Ndry_tot_iso = sum(Nisodry);
	xdryiso = bsxfun(@rdivide, Nisodry, Ndry_tot_iso);
	%xdryf = [Nif(:,1) Nif(:,2) Nif(:,3) Nif(:,5)] ./ Ntotf';

	figure 
	plot(vol * 10^3, x, '-r', 'LineWidth', 1.5), hold on
	plot(volu * 10^3, xdryiso, 'Color', 'blue', 'LineWidth', 1.5), hold on
	%plot(vola * 10^3, xdryf, '-g', 'LineWidth', 1.5), hold on
	title('Molar fractions profiles along the volume of catalyst');
	xlabel('Volume of catalyst [l]');
	ylabel('Molar fractions of species');

	%%
	% Define conversion plot
	conv_CH4 = (Ni0(3) - Ni(:, 3)) ./ (Ni0(3));
	conviso_CH4 = (Niso0(3) - Niso(:, 3)) ./ (Niso0(3));
	conv_CH4_f = (Ni0(3) - Nif(:, 3)) ./ (Ni0(3));

	figure
	plot(vol * 10^3, conv_CH4, '-r', 'LineWidth', 1.5), hold on
	plot(volu * 10^3, conviso_CH4, 'Color', 'blue', 'LineWidth', 1.5), hold on
	plot(vola * 10^3, conv_CH4_f, '-g', 'LineWidth', 1.5), hold on
	l = refline(0, 0.5);
	set(l, 'LineStyle', '--', 'Color', 'black');
	title('Conversion of methane along reactor volume');
	xlabel('Volume of catalyst[L]');
	ylabel('Conversion of Methane');

	%% 
	% Define yields plots
	yi_h2_ch4_iso = ((Niso(:, 2) - Niso0(2)) ./ Niso0(3)) .* (1/4);
	yi_h2_ch4 = ((Ni(:, 2) - Ni0(2)) ./ Ni0(3)) .* (1/4);
	yi_co_ch4_iso = (Niso(:, 1) ./ Niso0(3));
	yi_co_ch4 = (Ni(:, 1) ./ Ni0(3));
	yi_h2_ch4_f = ((Nif(:, 2) - Ni0(2)) ./ Ni0(3)) .* (1/4);
	yi_co_ch4_f = (Nif(:, 1) ./ Ni0(3));

	figure
	plot (volu * 10^3, yi_co_ch4_iso, 'Color', 'blue', 'LineWidth', 1.5), hold on, grid on
	plot(volu * 10^3, yi_h2_ch4_iso, 'Color', 'blue', 'LineWidth', 1.5)
	plot(vola * 10^3, yi_h2_ch4_f, '-g', 'LineWidth', 1.5)
	plot(vola * 10^3, yi_co_ch4_f, '-g', 'LineWidth', 1.5), hold on
	plot (vol * 10^3, yi_co_ch4, '-r', 'LineWidth', 1.5), hold on, grid on
	plot(vol * 10^3, yi_h2_ch4, '-r', 'LineWidth', 1.5), hold on
	xlabel('Volume of catalyst [l]');
	ylabel('Yield of H2 and CO on CH4');
	title('Yield of H2 and CO on CH4 along the volume of catalyst');


	%%
	% Define selectivities plot
	sh2 = (Ni(:, 2) - Ni0(2)) ./ (Ni0(3) - Ni(:, 3)) * 1/4;
	sco = (Ni(:, 1) - Ni0(1)) ./ (Ni0(3) - Ni(:, 3));
	sh2iso = (Niso(:, 2) - Niso0(2)) ./ (Niso0(3) - Niso(:, 3)) * 1/4;
	scoiso = (Niso(:, 1) - Niso0(1)) ./ (Niso0(3) - Niso(:, 3));
	scof = (Nif(:, 1) - Ni0(1)) ./ (Ni0(3) - Nif(:, 3));
	sh2f = (Nif(:, 2) - Ni0(2)) ./ (Ni0(3) - Nif(:, 3)) * 1/4;

	figure
	plot(vol * 10^3, sh2, '-r', 'LineWidth', 1.5),hold on, grid on
	plot(vol * 10^3, sco,'-r', 'LineWidth', 1.5) 
	plot(vola * 10^3, scof,'-g', 'LineWidth', 1.5) 
	plot(vola * 10^3, sh2f, '-g', 'LineWidth', 1.5), hold on, grid on
	plot(volu * 10^3, sh2iso, 'Color','blue', 'LineWidth', 1.5), hold on, grid on
	plot(volu * 10^3, scoiso,'Color','blue', 'LineWidth', 1.5), hold on
	title('Selectivities of H2,CO on CH4 along the volume of catalyst');
	xlabel('Volume of catalyst [l]');
	ylabel('Selectivity of H2,CO on CH4');

	%%
	% Residence time
	voliso = ((Niso * 3.6) .* (Rg ./ 1000 * T)) ./ P;
	volflow_totiso = voliso(:, 1) + voliso(:, 2) + voliso(:, 3) + voliso(:, 4) + voliso(:, 5);
	h = cumtrapz(volu, 1 ./ (volflow_totiso));

	volno = ((Ni .* 3.6) .* (Rg ./ 1000 .* Tk)) ./ P;
	volflow_tot = volno(:, 1) + volno(:, 2) + volno(:, 3) + volno(:, 4) + volno(:, 5);

	hno = cumtrapz(vol, 1 ./ (volflow_tot));

	volf = ((Nif .* 3.6) .* (Rg ./ 1000 .* Tkf)) ./ P;
	volflow_tot_f = volf(:, 1) + volf(:, 2) + volf(:, 3) + volf(:, 4) + volf(:, 5);

	hnof = cumtrapz(vola, 1 ./ (volflow_tot_f));

	figure
	plot(volu * 10^3, h * 10 ^ 3, 'Color', 'blue', 'LineWidth', 1.5), hold on
	plot(vol * 10^3, hno * 10 ^ 3, '-r', 'LineWidth', 1.5), hold on
	plot(vola * 10^3, hnof * 10 ^ 3, '-g', 'LineWidth', 1.5), hold on
	title('Residence time profile along the volume of catalyst');
	xlabel('Volume of catalyst [l]');
	ylabel('Residence time [ms]');

	%% 
	% Heat of reaction

	for u = 1:length(vol)
		[PFR, QR] = Adiabatic(vol(u), Y(u, 1:6)', nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai);
		Qvec(u) = QR;
	end

	for u = 1:length(volu)
		[PFR2, QR2] = Isothermal(volu(u), Y2(u, 1:6)', nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai);
		Qvec2(u) = QR2;
	end

	for u = 1:length(vola)
		[PFR3, QR3] = HeatFlux(vola(u), Y3(u, 1:6)', nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai);
		Qvec3(u) = QR3;
	end

	figure
	plot(vol * 10^3, Qvec, '-r', 'LineWidth', 1.5), hold on
	plot(volu * 10^3, Qvec2, 'Blue', 'LineWidth', 1.5), hold on
	plot(vola * 10^3, Qvec3, '-g', 'LineWidth', 1.5), hold on
	title('Heat of reaction vs volume of catalyst');
	xlabel('Volume of catalyst[l]');
	ylabel('Heat of reaction [KJ]');
	xlim([0 2]);
	ylim([0 15000]);

end
