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

%% Guess on parameters
par0 = [60 19500];

%% Variables declaration
global results_C1;
global results_t1;
Tvec_len = length(Tvec);
results_C1 = {};
results_t1 = {};
A_vec = zeros(1,Tvec_len);
Ea_vec = zeros(1,Tvec_len);
power_vec = zeros(1,Tvec_len);

%% Inlet and experimental data, search for the minimum of error function
for Tidx = 1 : Tvec_len
	
	temperature = Tvec(Tidx);
	time = Trials_trimmed(:,1);							% s
	coutsoda = Trials_trimmed(:, (Tidx + 1)) / 10^6;	% mol/cm^3;
	
	for i = 1 : length(coutsoda)
		if coutsoda(i) == 0
			time = time(1:(i - 1));
			coutsoda = coutsoda(1:(i - 1));
			break;
		end
	end
	cinsoda = coutsoda(1, 1);
	
	[par,fval] = fminsearch(@err, par0, opt, time, cinsoda, temperature, coutsoda, Tidx);
	
	A_vec(Tidx) = par(1);
	Ea_vec(Tidx) = par(2);
	
	a = plot(results_t1{Tidx},results_C1{Tidx},time,coutsoda,'o'),hold on
	set(a, 'Color', [0.2+(0.2*Tidx) 0.2 0.4], 'LineWidth', 1.25);
	title('Fitted experimental curves for a multiphase batch reactor under KC');
	ylabel('Concentration of OH- [mol/L]');
	xlabel('Time [s]');
	legend('calculated Cout_O_H_- at 5°C','Cexp,out_O_H_- at 5°C',...
		'calculated Cout_O_H_- at 17.5°C','Cexp,out_O_H_- at 17.5°C',...
		'calculated Cout_O_H_- at 30.6°C','Cexp,out_O_H_- at 30.6°C',...
		'calculated Cout_O_H_- at 43°C','Cexp,out_O_H_- at 43°C','Location','Northeast');
	drawnow
	
end

%% Fitted parameters for different temperatures
disp('Fitted pre-exponential factor: ');
pre_exp = mean(A_vec);
disp('Fitted activation energy: ');
act_energy = mean(Ea_vec);

%% Lovely message
disp('poop love');
	
%% Error function
function S = err(par, time, cinsoda, temperature,coutsoda, index)

	sol = Batch(par, time, cinsoda, temperature);
	Ccalc = deval(sol, time)';
	S = norm(Ccalc - coutsoda);
	t1   = linspace(min(time),max(time),100); 
	C1   = deval(sol,t1);
	global results_C1;
	global results_t1;
	results_t1{index} = t1;
	results_C1{index} = C1;
	
end

%% Solution function
function sol = Batch(par, time, cinsoda, temperature)

	nu = -1;
	Rg = 8.314;									% J/K mol
	porosity = 0.225;							% [-]
	soda = 470 * 10^3;							% cm^3
	resin = 10 * 10^3							% cm^3
	vol = soda + resin;							% cm^3
	xsoda = soda / vol;							% [-]
	xresin = resin / vol;						% [-]
	porosity = 0.225;							% [-]
	D = 650e-4;									% Resins sphere diameter - [cm]
	fiS = resin / vol;							% Solid fraction [-]
	fiL = soda / vol;							% Liquid fraction [-]
	aS = 6 / D;									% Solid specific area - [1/cm]
	aL = fiS * aS / fiL;						% Liquid specific area - [1/cm]
	
	sol = ode15s(@BMi, time, cinsoda, [], par, temperature, Rg, nu, aL);
	
end

%% Material balance
function Cprimo = BMi(time, c, par, Tin, Rg, nu, aL)

	A = par(1);
	Ea = par(2);
	k = A * exp(-Ea / Rg / Tin);
	R = k * c * aL;
	r = nu * R;
	Cprimo = r';
	
end













