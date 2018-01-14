clc,clear all,close all

% Reaction is R-plus + Naoh = R-na + oh-
%% Read experimental data from xls
Trials = xlsread('Expdata.xlsx');
Trials_trimmed = Trials(8:end, :);

%% constants
Rg = 8.314;										% J/K mol
nu = -1;

%% Data
soda = 470*1e-6;								% L
resin = 10*1e-6;								% L
vol = soda+resin;								% L
porosity = 0.225;								% [-]
D = 650e-4;										% Resins sphere diameter - [cm]
fiS = resin / vol;								% Solid fraction [-]
fiL = soda / vol;								% Liquid fraction [-]
aS = 6 / D;										% Specific solid area - [1/cm]
aL = fiS * aS / fiL;							% Specific liquid area - [1/cm]
Tvec = [ 5 17.5 30.6 43] + 273.15;				% K


%% Options
opt = optimset('Display','Iter');



%% Variables declaration
global results_C1;
global results_t1;
Tvec_len = length(Tvec);
results_C1 = {};
results_t1 = {};

%% Lovely message
disp('poop love');
	
%% Inlet data and experimental data check
for Tidx = 1 : Tvec_len
	
	temperature = Tvec(Tidx);
	time = Trials_trimmed(:,1);						% s
	coutsoda = Trials_trimmed(:, (Tidx + 1))/10^6;	% mol/cm^3;
	
		for i = 1 : length(coutsoda)
			if coutsoda(i) == 0
				time = time(1:(i - 1));
				coutsoda = coutsoda(1:(i - 1));
				break;
			end
		end
	cinsoda = coutsoda(1, 1);
	
	solution = Batch(time,cinsoda,temperature);
	
	a = plot(solution.x,solution.y,time,coutsoda,'o'),hold on
	set(a, 'Color', [1/(Tidx) 0 0], 'LineWidth', 1.25);
	title('Fitted experimental curves for a multiphase batch reactor under MTC');
	ylabel('Concentration of OH- [mol/L]');
	xlabel('Time [s]');
	legend('calculated Cout_O_H_- at 5°C','Cexp,out_O_H_- at 5°C',...
		'calculated Cout_O_H_- at 17.5°C','Cexp,out_O_H_- at 17.5°C',...
		'calculated Cout_O_H_- at 30.6°C','Cexp,out_O_H_- at 30.6°C',...
		'calculated Cout_O_H_- at 43°C','Cexp,out_O_H_- at 43°C','Location','Northeast');

end

%% Function solution
function sol = Batch(time, cinsoda, temperature)

	nu = -1;
	Rg = 8.314;														% J/K mol
	porosity = 0.225;												% [-]
	soda = 470 * 10^3;												% cm^3
	resin = 10 * 10^3												% cm^3
	vol = soda + resin;												% cm^3
	xsoda = soda / vol;												% [-]
	xresin = resin / vol;											% [-]
	porosity = 0.225;												% [-]
	D = 650e-4;														% Resins sphere diameter - [cm]
	fiS = resin / vol;												% Solid fraction [-]
	fiL = soda / vol;												% Liquid fraction [-]
	aS = 6 / D;														% Solid specific area - [1/cm]
	aL = fiS * aS / fiL;											% Liquid specific area - [1/cm]
	diff = Rg * 10^-3 * temperature /(96500 * (1/50.1 + 1/197.6));	% Diffusivity coefficient
	rho = 2.13 * 10^-3;												% Kg/cm^3
	vrel = 1;														% cm/s
	mu = 0.087														% Pa*s
	Rep = rho * vrel * D / mu;										% Reynolds number
	Sc = mu / rho / diff;											% Schmidt number
	Sh = 2 + 0.44 * Rep^0.5 * Sc^0.38;								% Sherwood number
	hm = Sh * diff / D;												% cm/s
	
	sol = ode15s(@BMi, time, cinsoda, [], temperature, Rg, nu, aL, hm);
	
end
	
function Cprimo = BMi(time, c, Tin, Rg, nu, aL, hm)

	R = hm * c * aL;
	r = nu * R;
	Cprimo = r';
	
end













