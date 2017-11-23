
% =========================== STEAM REFORMING =============================
%
% Calculates mass balances for a packed bed differential reactor.
% The molar flow rates are integrated upon the volume of catalyst.
% Plots molar flow rates profile along the reactor,
% dry products conversion, conversion of methane, selectivities,
% yields, gas velocities, residence time.
% Provides graphical indications of catalyst volume, yields, selectivities,
% gas velocities, residence time at 50% methane conversion.
%
% The listed reactions are a simplified version of methane steam reforming:
% 1) CH4 + H2O <=> CO + 3H2 (methane steam reforming)
% 2) CO + H2O <=> CO2 + H2 (water-gas shift)
% 3) CH4 + 2H2O <=> CO2 + 4H2 (methanation)
% In this application, the overall density is not constant and the reactor 
% is considered to be working at isothermal conditions.
%
% =========================================================================

function Isothermal_reactor

	%% ======================== MATERIAL BALANCE ==========================

	function [PFR] = BMi(vol, Ni, ki,kj, kpj, nu, phi, ro_bulk, P)
		% Define partial pressure as function of molar flow rates  
		Ntot = sum(Ni);
		p = P * Ni ./ Ntot;			% i = CO, H2, CH4, H2O, CO2

		% Calculate reaction rates
		DEN = 1 + ki(1) * p(1) + ki(2) * p(2) + ki(3) * p(3) + ((ki(4) * p(4)) / p(2)); 
		R1 = ((kj(1) / (p(2) ^ 2.5)) * (p(3) * p(4) - ((p(2) ^ 3 * p(1)) / kpj(1)))) / ((DEN) ^ 2);
		R2 = ((kj(2) / p(2)) * (p(1) * p(4) - ((p(2) * p(5)) / kpj(2)))) / ((DEN) ^ 2);
		R3 = ((kj(3) / (p(2) ^ 3.5)) * ((p(3) * (p(4) ^ 2)) - (p(2) ^ 4 * p(5) / kpj(3)))) / ((DEN) ^ 2);

		% Store reaction rates
		R = [R1 R2 R3]';

		% Calculate production rates
		r = nu * R;

		% Define material balance to integrate on the volume of catalyst
		PFR = r * (ro_bulk / phi);
	end


	%% =========================== MAIN SCRIPT ============================

	clc, clear all, close all;
	format long;

	% Initialize operating conditions
	T = 650 + 273.15;								% Kelvin
	P = 5;											% Bar

	% Initialize constants
	Rg = 8.3145 * 10^-2;							% Bar m3 / K kmol

	% Initialize kinetic and equilibrium parameters
	Aj = [4.225e+15 1.955e+06 1.020e+15]';			% With j = 1,2,3 - bar^-1
	Ai = [8.23e-05 6.12e-09 6.65e-04 1.77e+05]';	% With i = CO,H2,CH4,H2O - bar^-1, adim
	Ej = [240.1 67.13 243.9]' * 10;					% Bar m3/kmol
	DHi = [-70.65 -82.90 -38.28 88.68]' * 10;		% Bar m3/kmol
	ro_bulk = 1900;									% Kg/m3_reactor
	phi = 0.6;										% Adim

	% Calculate equilibrium constants, with j=1,2,3 - bar2, adim, bar2
	kpj = [exp(30.481 - 27.187e+03 / T) exp(-3.924 + 4.291e+03 / T) exp(26.891 - 23.258e+03 / T)]';

	% Calculate kinetic constants
	kj = Aj .* (exp(-Ej ./ (Rg * T)));
	ki = Ai .* (exp(-DHi ./ (Rg * T)));

	% Inlet molar flow rates
	Ni0 = ([0 2.63 50 150 0]' .* 3.6);				% From mol/s to kmol/h           

	% Define stoichiometric matrix
	nu=[+1 -1  0
		 3  1  4
		-1  0 -1
		-1 -1 -2
		 0  1  1];


	%% Integration
	[vol, Ni] = ode23(@BMi, [0 0.05], Ni0, [], ki, kj, kpj, nu, phi, ro_bulk, P);


	%% Plots

	% Define specie's flow rate profile along the volume of catalyst
	Nprof = [Ni(:, 1) Ni(:, 2) Ni(:, 3) Ni(:, 4) Ni(:, 5)];

	figure
	plot(vol * 10^3, Nprof, 'LineWidth', 1.5), hold on
	title('Flow rates profiles along the volume of catalyst');
	legend('CO', 'H2','CH4', 'H2O', 'CO2');
	xlabel('Volume of catalyst [l]');
	ylabel('Molar flow rate of species [kmol/h]');


	%%
	% Define dry products plot
	Ndry = [Ni(:, 1), Ni(:, 2), Ni(:, 3), Ni(:, 5)]';
	Ndry_tot = sum(Ndry);
	xdry = bsxfun(@rdivide, Ndry, Ndry_tot);

	figure
	plot(vol * 10^3, xdry, 'LineWidth', 1.5), hold on
	title('Dry products fraction along the volume of catalyst');
	legend('CO', 'H2', 'CH4', 'CO2');
	xlabel('Volume of catalyst [l]');
	ylabel('Molar fraction of dry products');


	%%
	% Define conversion plot
	conv_CH4 = (Ni0(3) - Ni(:, 3)) ./ (Ni0(3));

	figure
	plot(vol * 10^3, conv_CH4, 'LineWidth', 1.5), hold on, grid on
	l = refline(0, 0.5);
	set(l, 'LineStyle', '--', 'Color', 'black');
	l = line([4.18 4.18], [0 1]);
	set(l, 'LineStyle', '--', 'Color', 'black');
	title('Conversion of Methane along the volume of catalyst');
	xlabel('Volume of catalyst [l]');
	ylabel('Conversion of Methane');
	ylim([0 1]);


	%%
	% Define yields plots
	yi_h2_ch4 = ((Ni(:, 2) - Ni0(2)) ./ (Ni0(3))) .* (1/4);	% Defined on the last reaction
	yi_co_ch4 = (Ni(:, 1) ./ Ni0(3));

	figure
	plot(vol * 10^3, yi_co_ch4, 'LineWidth', 1.5), hold on, grid on
	plot(vol * 10^3, yi_h2_ch4, '-r', 'LineWidth', 1.5)
	l = line([4.18 4.18], [0 1]);
	set(l, 'LineStyle', '--', 'Color', 'black');
	title('Yield of H2 and CO on CH4 along the volume of catalyst');
	legend('CO', 'H2');
	xlabel('Volume of catalyst [l]');
	ylabel('Yield of H2 and CO on CH4');


	%%
	%Define selectivities plot
	sh2 = (Ni(:, 2) - Ni0(2)) ./ (Ni0(3) - Ni(:, 3)) * 1/4;	% Defined on the last reaction
	sco = (Ni(:, 1) - Ni0(1)) ./ (Ni0(3) - Ni(:, 3));

	figure
	plot(vol * 10^3, sh2, '-r', 'LineWidth', 1.5), hold on, grid on
	plot(vol * 10^3, sco, 'LineWidth', 1.5) 
	l = line([4.18 4.18], [0 1]);
	set(l, 'LineStyle', '--', 'Color', 'black');
	title('Selectivities of H2,CO on CH4 along the volume of catalyst');
	legend('H2', 'CO');
	xlabel('Volume of catalyst [l]');
	ylabel('Selectivity of H2,CO on CH4');


	%%
	% Ratio between inlet and outlet gas velocity
	vo = (Ni .* (Rg * T)) ./ P;
	volflow_out_tot = vo(:, 1) + vo(:, 2) + vo(:, 3) + vo(:, 4) + vo(:, 5);
	vi = (Ni0 .* (Rg * T)) ./ P;
	volflow_in_tot = vi(1) + vi(2) + vi(3) + vi(4) + vi(5);
	ratio = (volflow_out_tot / volflow_in_tot);

	figure
	plot(vol * 10^3, ratio, 'LineWidth', 1.5), hold on, grid on
	l = line([4.18 4.18], [1 1.5]);
	set(l, 'LineStyle', '--', 'Color', 'black');
	title('Outlet and inlet gas velocity ratio along the volume of catalyst');
	ylabel('Outlet gas velocity and inlet gas velocity ratio');
	xlabel('Volume of catalyst [l]');
	%ylim([1 1.35]);


	%%
	% Residence time plot
	volu = (Ni .* (Rg * T)) ./ P;
	volflow_tot = volu(:,1) + volu(:, 2) + volu(:, 3) + volu(:, 4) + volu(:, 5);
	h = cumtrapz(vol, 1 ./ (volflow_tot));

	figure
	plot(vol * 10^3, h * 3600 * 10^3, 'LineWidth', 1.5), grid on
	l = line([4.18 4.18], [0 4]);
	set(l, 'LineStyle', '--', 'Color', 'black');
	title('Residence time profile along the volume of catalyst');
	xlabel('Volume of catalyst [l]');
	ylabel('Residence time [ms]');
	ylim([0 2.5]);


	%%
	% Residence time at reactor outlet
	tau = trapz(vol, 1 ./ volflow_tot);

	format short
	disp('Residence time value at output [ms]: ');
	disp(tau * 3600 * 10^3);

end
