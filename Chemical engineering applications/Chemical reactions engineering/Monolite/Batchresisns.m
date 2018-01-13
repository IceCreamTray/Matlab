clc,clear all,close all
% Read experimental data
Trials = xlsread('Expdata.xlsx');
Trials_trimmed = Trials(8:end, :);

%constants
Rg = 8.314 %J/kmol
nu = -1;
% Reaction is R-plus + Naoh = R-na + oh-
% Data
t = Trials_trimmed(:,1);	%s
soda = 470*1e-3	%L
resin = 10*1e-3	%L
vol = soda+resin;	%L
xsoda = soda/vol;
xresin = resin/vol;
porosity = 0.225;	%[-]
Tin = 5 + 273.15	%K

% Concentrations
cinsoda = Trials_trimmed(1,2)/1e3;	%mol/L
coutsoda = Trials_trimmed(:,2)/1e3;
%plot to check
%scatter(t,coutsoda);

%axes color
set(0,'DefaultAxesColorOrder',[1 0 0; .1 1 .1; 0 0 1]);
opt=optimset('Display','Iter');
%parameters
par0 = [0.1 0.1];

[par,fval] = fminsearch(@err,par0,opt,t,coutsoda);   
disp(fval);
disp(par);

function S = err(par,t,coutsoda)
	sol = Batch(par,t);
	Ccalc = deval(sol,t)';
	S = norm(Ccalc - coutsoda);
	t1   = linspace(min(t),max(t),100); 
	C1   = deval(sol,t1);
	plot(t1,C1,t,coutsoda,'s')
	drawnow
end

function sol = Batch(par,t)
	Trials = xlsread('Expdata.xlsx');
	Trials_trimmed = Trials(8:end, :);
	cinsoda = Trials_trimmed(1,2)/1e3;
	nu = -1;
	Tin = 5 + 273.15;	%K
	Rg = 8.314; %J/kmol
	sol = ode15s(@BMi,t,cinsoda,[],par,Tin,Rg,nu);
end
function Cprimo = BMi(t,c,par,Tin,Rg,nu)
	A=par(1);
	Ea=par(2);
	k = A*exp(-Ea/Rg/Tin);
	R = k*c;
	r = nu*R;
	Cprimo = r';
end













