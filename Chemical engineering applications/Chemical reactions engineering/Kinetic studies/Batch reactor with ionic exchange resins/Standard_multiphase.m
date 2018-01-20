%clc,clear all,close all

% Reaction is R-plus + Naoh = R-na + oh-
%% Read experimental data from xls
Trials = xlsread('Expdata.xlsx');
Trials_trimmed = Trials(8:end, :);

%% constants
Rg = 8.314;														% J/K mol
nu = -1;

%% Data
Tvec = [ 5 17.5 30.6 43] + 273.15;

%% Colors
Lcol = { [1 0 0] [0 1 0] [0 0 1] [0 0 0] };

%% Guess on parameters
par0 = [50 10000];

%% Variables declaration
global results_C1;
global results_t1;
global hm;
Tvec_len = length(Tvec);
results_C1 = {};
results_t1 = {};
A_vec = zeros(1,Tvec_len);
Ea_vec = zeros(1,Tvec_len);

%% Options
opt = optimset('Display','Iter');


%% Inlet and experimental data, search for the minimum of error function
for Tidx = 1 : Tvec_len
	
	temperature = Tvec(Tidx);
	time = Trials_trimmed(:,1);									% s
	coutsoda = Trials_trimmed(:, (Tidx + 1)) / 10^6;			% mol/cm^3;
	
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

	a = plot(results_t1{Tidx},results_C1{Tidx},time,coutsoda,'o');hold on
	set(a, 'Color', Lcol{Tidx}, 'LineWidth', 1.25);
	title('Fitted experimental curves for a multiphase batch reactor, general case');
	ylabel('Concentration of OH- [mol/L]');
	xlabel('Time [s]');
	legend('calculated Cout_O_H_- at 5�C','Cexp,out_O_H_- at 5�C',...
		'calculated Cout_O_H_- at 17.5�C','Cexp,out_O_H_- at 17.5�C',...
		'calculated Cout_O_H_- at 30.6�C','Cexp,out_O_H_- at 30.6�C',...
		'calculated Cout_O_H_- at 43�C','Cexp,out_O_H_- at 43�C','Location','Northeast');
	drawnow
	
end

%% Fitted parameters for different temperatures
disp('Fitted pre-exponential factor: ');
pre_exp = mean(A_vec)
disp('Fitted activation energy: ');
act_energy = mean(Ea_vec)

%% Different k evaluation
figure
k1 = log(A_vec(1)*exp(-Ea_vec(1)/Rg/Tvec(1)));
k2 = log(A_vec(2)*exp(-Ea_vec(2)/Rg/Tvec(2)));
k3 = log(A_vec(3)*exp(-Ea_vec(3)/Rg/Tvec(3)));
k4 = log(A_vec(4)*exp(-Ea_vec(4)/Rg/Tvec(4)));
scatter(1./Tvec,[k1 k2 k3 k4]),hold on;
p = polyfit(1./Tvec,[k1 k2 k3 k4],1);
f = polyval(p,1./Tvec);
plot(1./Tvec,f,'Color','Black');
title('Control on calculated kinetic constants');
ylabel('kinetic constants');
xlabel('1/T [1/K]');

%% Arrhenius plot
figure
k = pre_exp*exp(-act_energy/Rg./Tvec);
kapp = k.*hm./(hm+k);
lnk = log(kapp);
plot(1./Tvec,lnk, 'Linewidth', 1.5);
title('Kinetic constant as function of T');
ylabel('lnK');
xlabel('1/T [1/K]');

%% Lovely message
%disp('poop love');

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

function sol = Batch(par, time, cinsoda, temperature)

	nu = -1;
	Rg = 8.314;														% J/K mol
	porosity = 0.225;												% [-]
	soda = 470 * 10^3;												% cm^3
	resin = 10 * 10^3;												% cm^3
	vol = soda + resin;												% cm^3
	xsoda = soda / vol;												% [-]
	xresin = resin / vol;											% [-]
	porosity = 0.225;												% [-]
	D = 650e-4;														% Resins sphere diameter - [cm]
	fiS = resin /vol*porosity;										% Solid fraction [-]
	fiL = soda / vol;												% Liquid fraction [-]
	aS = 6 / D;														% Solid specific area - [1/cm]
	aL = fiS * aS / fiL;											% Liquid specific area - [1/cm]
	diff = Rg * 10^-3 * temperature /(96500 *...
			(1/50.1 + 1/197.6))*10^3;								% Diffusivity coefficient
	rho = 2.13 * 10^-3;												% Kg/cm^3
	vrel = 1;														% cm/s
	mu = 0.087;														% Pa*s
	Rep = rho * vrel * D / mu;										% Reynolds number
	Sc = mu / rho / diff;											% Schmidt number
	Sh = 2 + 0.44 * Rep^0.5 * Sc^0.38;								% Sherwood number
	global hm
	hm = Sh * diff / D;												% cm/s
	
	sol = ode15s(@BMi, time, cinsoda, [], par, temperature, Rg, nu, aL, hm);
	
end
	
function Cprimo = BMi(time,c,par,Tin,Rg,nu,aL,hm)

	A = par(1);
	Ea = par(2);
	k = A*exp(-Ea/Rg/Tin);
	kapp = hm*k/(hm+k);
	R = kapp*c*aL;
	r = nu*R;
	Cprimo = r';
	
end






