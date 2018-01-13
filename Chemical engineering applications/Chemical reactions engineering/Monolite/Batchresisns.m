clc,clear all,close all
% Read experimental data
Trials = xlsread('Expdata.xlsx');
Trials_trimmed = Trials(8:end, :);

%constants
Rg = 8.314 %J/kmol
nu = -1;
% Reaction is R-plus + Naoh = R-na + oh-
% Data
soda = 470*1e-3	%L
resin = 10*1e-3	%L
vol = soda+resin;	%L
xsoda = soda/vol;
xresin = resin/vol;
porosity = 0.225;	%[-]

% Concentrations
%plot to check
%scatter(t,coutsoda);

%axes color
set(0,'DefaultAxesColorOrder',[1 0 0; .1 1 .1; 0 0 1]);
opt=optimset('Display','Iter');

par0 = [0.1 0.1];


global results_C1;
global results_t1;
results_C1 = {};
results_t1 = {};

Tvec = [ 5 17.5 30.6 43] + 273.15;
Tvec_len = length(Tvec);

time = Trials_trimmed(:,1);	%s


[par,fval] = fminsearch(@err, par0, opt, Trials_trimmed, Tvec, Tvec_len, time);
disp(fval);
disp(par);

function S = err(par, data, Tvec, Tvec_len, time)
	global results_C1;
	global results_t1;
	
	for Tidx = 1 : Tvec_len
		temperature = Tvec(Tidx);
		coutsoda = data(:, (Tidx + 1)) / 1e3;
		cinsoda = coutsoda(1, 1);
	
		sol = Batch(par, time, cinsoda, temperature);
		Ccalc = deval(sol, time)';
		S = norm(Ccalc - coutsoda);
		t1   = linspace(min(time),max(time),100); 
		C1   = deval(sol,t1);
		
		results_t1{Tidx} = t1;
		results_C1{Tidx} = C1;
	end
end

for i = 1 : Tvec_len
	plot(results_t1{i}, results_C1{i}, time, Trials_trimmed(:, (i + 1)), 's'), hold on
	drawnow
end

function sol = Batch(par, time, cinsoda, temperature)
	nu = -1;
	Rg = 8.314; %J/kmol
	sol = ode15s(@BMi,time,cinsoda,[],par,temperature,Rg,nu);
end
	
function Cprimo = BMi(time,c,par,Tin,Rg,nu)
	A=par(1);
	Ea=par(2);
	k = A*exp(-Ea/Rg/Tin);
	R = k*c;
	r = nu*R;
	Cprimo = r';
end













