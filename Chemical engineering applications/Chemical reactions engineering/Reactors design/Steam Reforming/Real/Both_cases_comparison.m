% =========================== STEAM REFORMING =============================
%
% Calculates mass and energy balances for a packed bed real reactor.
% The molar flow rates are integrated upon the volume of catalyst.
% Plots molar flow rates profile along the reactor, conversions, yields,
% selectivities, heat of reaction, residence time, reaction rates.
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
	
	function[PFR, QR, Rd, Ri] = Adiabatic(vol, Y, nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai)
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
        
        Rd1 = ((kj(1) / (p(2) ^ 2.5)) * p(3) * p(4)) / (DEN ^ 2);
        Ri1 = ((kj(1) / (p(2) ^ 2.5)) * p(2) ^3 * p(1) / kpj(1)) / (DEN ^ 2);
        Rd2 = ((kj(2) / p(2)) * p(1) * p(4)) / (DEN ^ 2);
        Ri2 = ((kj(2) / p(2)) * p(2) * p(5) / kpj(2)) / (DEN ^ 2);
        Rd3 = ((kj(3) / p(2) ^ 3.5) * p(3) * p(4) ^ 2) / (DEN ^ 2);
        Ri3 = ((kj(3) / p(2) ^ 3.5) * p(2) ^4 * p(5) / kpj(3)) / (DEN ^ 2);
        
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
        
        Rd = [Rd1 Rd2 Rd3]';
        Ri = [Ri1 Ri2 Ri3]';

	end

	%% ========================= ISOTHERMAL CASE ==========================

	function[PFR2, QR2, Rdk2, Rik2] = Isothermal(volu, Y, nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai)
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
        
        Rd12 = ((kj(1) / (p(2) ^ 2.5)) * p(3) * p(4)) / (DEN ^ 2);
        Ri12 = ((kj(1) / (p(2) ^ 2.5)) * p(2) ^3 * p(1) / kpj(1)) / (DEN ^ 2);
        Rd22 = ((kj(2) / p(2)) * p(1) * p(4)) / (DEN ^ 2);
        Ri22 = ((kj(2) / p(2)) * p(2) * p(5) / kpj(2)) / (DEN ^ 2);
        Rd32 = ((kj(3) / p(2) ^ 3.5) * p(3) * p(4) ^ 2) / (DEN ^ 2);
        Ri32 = ((kj(3) / p(2) ^ 3.5) * p(2) ^4 * p(5) / kpj(3)) / (DEN ^ 2);
        
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
        
        Rdk2 = [Rd12 Rd22 Rd32]';
        Rik2 = [Ri12 Ri22 Ri32]';

	end

	%% ========================= HEAT FLOW CONST ==========================
	
	function [PFR3, QR3, Rda, Ria] = HeatFlux(vola, Y, nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai)
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
        
        Rd12 = ((kj(1) / (p(2) ^ 2.5)) * p(3) * p(4)) / (DEN ^ 2);
        Ri12 = ((kj(1) / (p(2) ^ 2.5)) * p(2) ^3 * p(1) / kpj(1)) / (DEN ^ 2);
        Rd22 = ((kj(2) / p(2)) * p(1) * p(4)) / (DEN ^ 2);
        Ri22 = ((kj(2) / p(2)) * p(2) * p(5) / kpj(2)) / (DEN ^ 2);
        Rd32 = ((kj(3) / p(2) ^ 3.5) * p(3) * p(4) ^ 2) / (DEN ^ 2);
        Ri32 = ((kj(3) / p(2) ^ 3.5) * p(2) ^4 * p(5) / kpj(3)) / (DEN ^ 2);
        
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
        
        Rda = [Rd12 Rd22 Rd32]';
        Ria = [Ri12 Ri22 Ri32]';

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
	blue = [0 128/256 255/256];
	red = [255/256 64/256 0];
	green = [128/256 224/256 128/256];

	%%
	% Profile along the volume of catalyst
	Nprof = [Ni(:, 1) Ni(:, 2) Ni(:, 3) Ni(:, 4) Ni(:, 5)];
	Nprofiso = [Niso(:, 1) Niso(:, 2) Niso(:, 3) Niso(:, 4) Niso(:, 5)];
	%Nproff = [Nif(:,1) Nif(:,2) Nif(:,3) Nif(:,4) Nif(:,5)];

	figure
	plot(vol * 10^3, Nprof, '-o', 'Color', red, 'LineWidth', 1.5), hold on
	plot(volu * 10^3, Nprofiso, '-o', 'Color',blue, 'LineWidth', 1.5), hold on
	%plot(vola * 10^3, Nproff, '-g', 'LineWidth', 1.5), hold on
	title('Flow rate profile along the volume of catalyst');
	xlabel('Volume of catalyst [l]');
	ylabel('Molar flow rate of species [kmol/h]');
	xlim([0 20]);

	%% 
	% Define dry molar fractions plot
	x = bsxfun(@rdivide, [Ni(:, 1) Ni(:, 2) Ni(:, 3) Ni(:, 5)], Nitot');
	Nisodry = [Niso(:, 1), Niso(:, 2), Niso(:, 3), Niso(:, 5)]';
	Ndry_tot_iso = sum(Nisodry);
	xdryiso = bsxfun(@rdivide, Nisodry, Ndry_tot_iso);
	%xdryf = [Nif(:,1) Nif(:,2) Nif(:,3) Nif(:,5)] ./ Ntotf';

	figure 
	plot(vol * 10^3, x, '-o', 'Color', red, 'LineWidth', 1.5), hold on
	plot(volu * 10^3, xdryiso,'-o', 'Color', blue, 'LineWidth', 1.5), hold on
	%plot(vola * 10^3, xdryf, '-g', 'LineWidth', 1.5), hold on
	title('Molar fractions profiles along the volume of catalyst');
	xlabel('Volume of catalyst [l]');
	ylabel('Molar fractions of species');
	xlim([0 20]);
	
	%%
	% Define temperature plot
	
	figure
	plot(vol * 10^3, Tk - 273.15, 'Color', red, 'LineWidth', 1.5), hold on
	o = refline(0, 650);
	set(o, 'Color', blue, 'LineWidth', 1.5);
	plot(vola * 10^3, Tkf - 273.15, 'Color', green, 'LineWidth', 1.5), hold on
	title('Temperature profiles along the volume of catalyst');
	xlabel('Volume of catalyst [l]');
	ylabel('Temperature [°C]');
	xlim([0 30])
	
	%%
	% Define conversion plot
	conv_CH4 = (Ni0(3) - Ni(:, 3)) ./ (Ni0(3));
	conviso_CH4 = (Niso0(3) - Niso(:, 3)) ./ (Niso0(3));
	conv_CH4_f = (Ni0(3) - Nif(:, 3)) ./ (Ni0(3));

	figure
	plot(vol * 10^3, conv_CH4, 'Color', red, 'LineWidth', 1.5), hold on
	plot(volu * 10^3, conviso_CH4, 'Color', blue, 'LineWidth', 1.5), hold on
	plot(vola * 10^3, conv_CH4_f, 'Color', green, 'LineWidth', 1.5), hold on
	l = refline(0, 0.5);
	set(l, 'LineStyle', '--', 'Color', 'black');
	title('Conversion of methane along reactor volume');
	xlabel('Volume of catalyst[L]');
	ylabel('Conversion of Methane');
	xlim([0 30]);

	%% 
	% Define yields plots
% 	yi_h2_ch4_iso = ((Niso(:, 2) - Niso0(2)) ./ Niso0(3)) .* (1/4);
% 	yi_h2_ch4 = ((Ni(:, 2) - Ni0(2)) ./ Ni0(3)) .* (1/4);
% 	yi_co_ch4_iso = (Niso(:, 1) ./ Niso0(3));
% 	yi_co_ch4 = (Ni(:, 1) ./ Ni0(3));
% 	yi_h2_ch4_f = ((Nif(:, 2) - Ni0(2)) ./ Ni0(3)) .* (1/4);
% 	yi_co_ch4_f = (Nif(:, 1) ./ Ni0(3));
% 
% 	figure
% 	plot (volu * 10^3, yi_co_ch4_iso, 'Color', blue, 'LineWidth', 1.5), hold on, grid on
% 	plot(volu * 10^3, yi_h2_ch4_iso, 'Color', blue, 'LineWidth', 1.5)
% 	plot(vola * 10^3, yi_h2_ch4_f, 'Color', green, 'LineWidth', 1.5)
% 	plot(vola * 10^3, yi_co_ch4_f, 'Color', green, 'LineWidth', 1.5), hold on
% 	plot (vol * 10^3, yi_co_ch4, 'Color', red, 'LineWidth', 1.5), hold on, grid on
% 	plot(vol * 10^3, yi_h2_ch4, 'Color', red, 'LineWidth', 1.5), hold on
% 	xlabel('Volume of catalyst [l]');
% 	ylabel('Yield of H2 and CO on CH4');
% 	title('Yield of H2 and CO on CH4 along the volume of catalyst');


	%%
	% Define selectivities plot
% 	sh2 = (Ni(:, 2) - Ni0(2)) ./ (Ni0(3) - Ni(:, 3)) * 1/4;
% 	sco = (Ni(:, 1) - Ni0(1)) ./ (Ni0(3) - Ni(:, 3));
% 	sh2iso = (Niso(:, 2) - Niso0(2)) ./ (Niso0(3) - Niso(:, 3)) * 1/4;
% 	scoiso = (Niso(:, 1) - Niso0(1)) ./ (Niso0(3) - Niso(:, 3));
% 	scof = (Nif(:, 1) - Ni0(1)) ./ (Ni0(3) - Nif(:, 3));
% 	sh2f = (Nif(:, 2) - Ni0(2)) ./ (Ni0(3) - Nif(:, 3)) * 1/4;
% 
% 	figure
% 	plot(vol * 10^3, sh2, 'Color', red, 'LineWidth', 1.5),hold on, grid on
% 	plot(vol * 10^3, sco,'Color', red, 'LineWidth', 1.5) 
% 	plot(vola * 10^3, scof,'Color', green, 'LineWidth', 1.5) 
% 	plot(vola * 10^3, sh2f, 'Color', green, 'LineWidth', 1.5), hold on, grid on
% 	plot(volu * 10^3, sh2iso, 'Color', blue, 'LineWidth', 1.5), hold on, grid on
% 	plot(volu * 10^3, scoiso,'Color', blue, 'LineWidth', 1.5), hold on
% 	title('Selectivities of H2,CO on CH4 along the volume of catalyst');
% 	xlabel('Volume of catalyst [l]');
% 	ylabel('Selectivity of H2,CO on CH4');

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
	plot(volu * 10^3, h * 10 ^ 3, 'Color', blue, 'LineWidth', 1.5), hold on
	plot(vol * 10^3, hno * 10 ^ 3, 'Color', red, 'LineWidth', 1.5), hold on
	plot(vola * 10^3, hnof * 10 ^ 3, 'Color', green, 'LineWidth', 1.5), hold on
	title('Residence time profile along the volume of catalyst');
	xlabel('Volume of catalyst [l]');
	ylabel('Residence time [ms]');
	xlim([0 30]);

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
	plot(vol * 10^3, Qvec, 'Color', red, 'LineWidth', 1.5), hold on
	plot(volu * 10^3, Qvec2, 'Color', blue, 'LineWidth', 1.5), hold on
	plot(vola * 10^3, Qvec3, 'Color', green, 'LineWidth', 1.5), hold on
	title('Heat of reaction vs volume of catalyst');
	xlabel('Volume of catalyst [l]');
	ylabel('Heat of reaction [kW]');
	xlim([0 30]);
	ylim([0 10000]);
    
    %%
    % Reaction rates
    for u = 1:length(vol)
		[PFR, QR, Rd, Ri] = Adiabatic(vol(u), Y(u, 1:6)', nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai);
        Rdir(:,u) = Rd;
        Rinv(:,u) = Ri;
    end
    Rdir1 = Rdir(1,1:101);
    Rdir2 = Rdir(2,1:101);
    Rdir3 = Rdir(3,1:101);
    Rinv1 = Rinv(1,1:101);
    Rinv2 = Rinv(2,1:101);
    Rinv3 = Rinv(3,1:101);
    
    figure
    plot(vol * 10^3, Rdir1, 'LineWidth', 1.5), hold on, grid on 
    plot(vol * 10^3, Rinv1, 'LineWidth', 1.5);
    title('Reaction rate of methane steam reforming [Adiabatic]');
    xlabel('Volume of catalyst [l]');
    ylabel('Reaction rate [kmol/Kg_cat h]');
    legend('Direct reaction','Inverse reaction');
    ylim([0 50]);
    xlim([0 20]);
    
    figure
    plot(vol * 10^3, Rdir2, 'LineWidth', 1.5), hold on, grid on 
    plot(vol * 10^3, Rinv2, 'LineWidth', 1.5);
    title('Reaction rate of water - gas shift [Adiabatic]');
    xlabel('Volume of catalyst [l]');
    ylabel('Reaction rate [kmol/Kg_cat h]');
    legend('Direct reaction','Inverse reaction');
    ylim([0 50]);
    xlim([0 20]);
    
    figure
    plot(vol * 10^3, Rdir3, 'LineWidth', 1.5), hold on, grid on 
    plot(vol * 10^3, Rinv3, 'LineWidth', 1.5);
    title('Reaction rate of methanation [Adiabatic]');
    xlabel('Volume of catalyst [l]');
    ylabel('Reaction rate [kmol/Kg_cat h]');
    legend('Direct reaction','Inverse reaction');
    ylim([0 50]);
    xlim([0 20]);
    
    for u = 1:length(volu)
		[PFR, QR, Rdk2, Rik2] = Isothermal(volu(u), Y2(u, 1:6)', nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai);
         Rdirk2(:,u) = Rdk2;
         Rinvk2(:,u) = Rik2;
    end
    Rdir21 = Rdirk2(1,1:94);
    Rdir22 = Rdirk2(2,1:94);
    Rdir23 = Rdirk2(3,1:94);
    Rinv21 = Rinvk2(1,1:94);
    Rinv22 = Rinvk2(2,1:94);
    Rinv23 = Rinvk2(3,1:94);
    rnet=Rdir21-Rinv21;
    figure
	plot(volu*10^3,rnet),hold on
    plot(volu * 10^3, Rdir21, 'LineWidth', 2), hold on, grid on 
    plot(volu * 10^3, Rinv21, 'LineWidth', 2);
    title('Reaction rate of methane steam reforming [Isothermal]');
    xlabel('Volume of catalyst [l]');
    ylabel('Reaction rate [kmol/Kg_c h]');
    legend('Direct reaction','Inverse reaction');
    ylim([0 50]);
    xlim([0 20]);
    
    figure
    plot(volu * 10^3, Rdir22, 'LineWidth', 2), hold on, grid on 
    plot(volu * 10^3, Rinv22, 'LineWidth', 2);
    title('Reaction rate of water - gas shift [Isothermal]');
    xlabel('Volume of catalyst [l]');
    ylabel('Reaction rate [kmol/Kg_c h]');
    legend('Direct reaction','Inverse reaction');
    ylim([0 50]);
    xlim([0 20]);
    
    figure
    plot(volu * 10^3, Rdir23, 'LineWidth', 2), hold on, grid on 
    plot(volu * 10^3, Rinv23, 'LineWidth', 2);
    title('Reaction rate of methanation [Isothermal]');
    xlabel('Volume of catalyst [l]');
    ylabel('Reaction rate [kmol/Kg_c h]');
    legend('Direct reaction','Inverse reaction');
    ylim([0 50]);
    xlim([0 20]);
    
     for u = 1:length(vola)
		[PFR, QR, Rda, Ria] = HeatFlux(vola(u), Y3(u, 1:6)', nu, phi, ro_bulk, P, Rg, ai, bi, DHi, Ej, Aj, Ai);
         Rdira(:,u) = Rda;
         Rinva(:,u) = Ria;
    end
    Rdira1 = Rdira(1,1:106);
    Rdira2 = Rdira(2,1:106);
    Rdira3 = Rdira(3,1:106);
    Rinva1 = Rinva(1,1:106);
    Rinva2 = Rinva(2,1:106);
    Rinva3 = Rinva(3,1:106);
    
    figure
    plot(vola * 10^3, Rdira1, 'LineWidth', 1.5), hold on, grid on 
    plot(vola * 10^3, Rinva1, 'LineWidth', 1.5);
    title('Reaction rate of methane steam reforming [Constant heat flow]');
    xlabel('Volume of catalyst [l]');
    ylabel('Reaction rate [kmol/Kg_cat h]');
    legend('Direct reaction','Inverse reaction');
    ylim([0 30]);
    xlim([0 20]);
    
    figure
    plot(vola * 10^3, Rdira2, 'LineWidth', 1.5), hold on, grid on 
    plot(vola * 10^3, Rinva2, 'LineWidth', 1.5);
    title('Reaction rate of water - gas shift [Constant heat flow]');
    xlabel('Volume of catalyst [l]');
    ylabel('Reaction rate [kmol/Kg_cat h]');
    legend('Direct reaction','Inverse reaction');
    ylim([0 30]);
    xlim([0 20]);
    
    figure
    plot(vola * 10^3, Rdira3, 'LineWidth', 1.5), hold on, grid on 
    plot(vola * 10^3, Rinva3, 'LineWidth', 1.5);
    title('Reaction rate of methanation [Constant heat flow]');
    xlabel('Volume of catalyst [l]');
    ylabel('Reaction rate [kmol/Kg_cat h]');
    legend('Direct reaction','Inverse reaction');
    ylim([0 30]);
    xlim([0 20]);

end
