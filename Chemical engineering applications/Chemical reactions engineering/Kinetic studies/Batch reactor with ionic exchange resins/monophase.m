clc,clear all,close all
% Reaction is R-plus + Naoh = R-na + oh-

%% Read experimental data from xlsx
Trials = xlsread('Expdata.xlsx');
Trials_trimmed = Trials(8:end, :);

%% constants
Rg = 8.314;										% J/K mol
nu = -1;							

%% Data
soda = 470 * 1e-3;									% L
resin = 10 * 1e-3;									% L
vol = soda + resin;									% L
phiL = soda / vol;									% [-]
phiS = resin / vol;									% [-]
porosity = 0.225;									% [-]
Tvec = [5 17.5 30.6 43] + 273.15;					% [K]

%% Colors
Lcol = { [1 0 0] [0 1 0] [0 0 1] [0 0 0] };

%% options
opt = optimset('Display','Iter');

%% Guess on parameters
par0 = [60 19500];

%% Variables definition
global results_C1;
global results_t1;
Tvec_len = length(Tvec);
results_C1 = {};
results_t1 = {};
A_vec = zeros(1,Tvec_len);
Ea_vec = zeros(1,Tvec_len);

%% Inlet and experimental data, search for the minimum of error function
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
	
	[par,fval] = fminsearch(@err, par0, opt, time, cinsoda, temperature, coutsoda, Tidx);
	
	A_vec(Tidx) = par(1);
	Ea_vec(Tidx) = par(2);
	

	a = plot(results_t1{Tidx},results_C1{Tidx}/10^-3,time,coutsoda/10^-3,'o');hold on
	set(a, 'Color', Lcol{Tidx}, 'LineWidth', 1.25);
	title('Fitted experimental curves for a monophase batch reactor')
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
pre_exp = mean(A_vec)
disp('Fitted activation energy: ');
act_energy = mean(Ea_vec)

%% Arrhenius plot
figure
k = pre_exp*exp(-act_energy/Rg./Tvec);
lnk = log(k);
plot(1./Tvec,lnk, 'Linewidth', 1.5);
title('Kinetic constant as function of T');
ylabel('lnK');
xlabel('1/T [1/K]');

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
	Rg = 8.314;													% J/ K mol
	sol = ode15s(@BMi,time,cinsoda,[],par,temperature,Rg,nu);
	
end

%% Mass balance for a monophase batch reactor, first order kinetics
function Cprimo = BMi(time,c,par,Tin,Rg,nu)

	A = par(1);
	Ea = par(2);
	k = A*exp(-Ea/Rg/Tin);
	R = k*c;
	r = nu*R;
	Cprimo = r';
	
end













