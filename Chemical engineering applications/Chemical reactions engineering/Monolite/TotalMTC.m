clc,clear all,close all
% Read experimental data
Trials = xlsread('Expdata.xlsx');
Trials_trimmed = Trials(8:end, :);

%constants
Rg = 8.314 %J/K mol
nu = -1;
% Reaction is R-plus + Naoh = R-na + oh-
% Data
soda = 470*1e-6	%m^3
resin = 10*1e-6	%m^3
vol = soda+resin;	%m^3
porosity = 0.225;	%[-]
D = 650e-6;             % Diametro delle resine - [m]
fiS = resin/vol;        % Frazione di solido [-]
fiL = soda/vol;        % Frazione di liquido [-]
aS = 6/D;               % Area specifica del solido - [1/m]
aL = fiS*aS/fiL;        % Area specifica del liquido - [1/m]

% Concentrations
%plot to check
%scatter(t,coutsoda);

%axes color
set(0,'DefaultAxesColorOrder',[1 0 0; .1 1 .1; 0 0 1]);
opt=optimset('Display','Iter');

par0 = 1;


global results_C1;
global results_t1;

Tvec = [ 5 17.5 30.6 43] + 273.15;
Tvec_len = length(Tvec);
results_C1 = {};
results_t1 = {};
A_vec = zeros(1,Tvec_len);
Ea_vec = zeros(1,Tvec_len);
power_vec = zeros(1,Tvec_len);
for Tidx = 1 : Tvec_len
	temperature = Tvec(Tidx);
	time = Trials_trimmed(:,1);	%s
	coutsoda = Trials_trimmed(:, (Tidx + 1))/10^6 %mol/cm^3;
	for i = 1 : length(coutsoda)
		if coutsoda(i) == 0
			time = time(1:(i - 1));
			coutsoda = coutsoda(1:(i - 1));
			break;
		end
	end
	cinsoda = coutsoda(1, 1);
	
	[par,fval] = fminsearch(@err, par0, opt, time, cinsoda, temperature, coutsoda, Tidx);
	disp(par)
	%power_vec(Tidx) = par(3);
	plot(results_t1{Tidx},results_C1{Tidx},time,coutsoda,'s'),hold on
	drawnow
end



disp('poop love');
	

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
% for Tidx = 1 : Tvec_len
% temperature = Tvec(Tidx);
% time = Trials_trimmed(:,1);	%s
% coutsoda = Trials_trimmed(:, (Tidx + 1))/10^6; %mol/cm^3;
% 	for i = 1 : length(coutsoda)
% 		if coutsoda(i) == 0
% 			time = time(1:(i - 1));
% 			coutsoda = coutsoda(1:(i - 1));
% 			break;
% 		end
% 	end
% cinsoda = coutsoda(1, 1);
% solution = Batch(time,cinsoda,temperature);
% plot(solution.x,solution.y),hold on
% scatter(time,coutsoda)
% end
function sol = Batch(par,time, cinsoda, temperature)
	nu = -1;
	Rg = 8.314; %J/ K mol
	porosity = 0.225;	%[-]
	soda = 470*10^3	%cm^3
	resin = 10*10^3	%cm^3
	vol = soda+resin;	%cm^3
	xsoda = soda/vol;
	xresin = resin/vol;
	porosity = 0.225;	%[-]
	D = 650e-4;             % Diametro delle resine - [cm]
	aS = 6/D;               % Area specifica del solido - [1/cm]
	aL = aS*(1-porosity)/porosity;
	mu = 0.087;	%cpoise
	diff = Rg*10^-3*temperature/(96500*(1/50.1+1/197.6));
	rho = 2.13*10^-3; %Kg/cm^3
	vrel = 1; %cm/s
	Rep = rho * vrel * D / mu;
	Sc = mu/rho/diff;
	Sh = 2 + 0.44*Rep^0.5*Sc^0.38;
	hm = Sh * diff/D; %cm/s
	sol = ode15s(@BMi,time,cinsoda,[],par,temperature,Rg,nu,aL,hm);
end
	
function Cprimo = BMi(time,c,par,Tin,Rg,nu,aL,hm)
	%k = A*exp(-Ea/Rg/Tin);
	R = hm*c;
	r = nu*R;
	Cprimo = r';
end













