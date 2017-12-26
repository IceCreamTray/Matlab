% Multiphase balances
function aaaa
clc,clear all
close all

% dati da gas
D   = 1e-5;   % m^2/s
v   = 1;      % m/s
rho = 1;    % kg/m3
mu  = 1e-5; % Pa*s

L = 8.5*10^-3;      % lunghezza reattore, m
Sc  =  mu/(rho*D);

d = 8e-3;       % m, diametro canale
a = 4/d;        % 1/m  area specifica
        
        %calcolo coeff di mass transfer 
Sh  =  3;
hm  =  Sh*D/d;  %m/s

M = [1 0
     0 0]; 
 
 % Read the experimental data from xls
%sheet1 = xlsread('Expdata.xlsx',1);
%sheet2 = xlsread('Expdata.xlsx',2);
%sheet3 = xlsread('Expdata.xlsx',3);
sheet4 = xlsread('Expdata.xlsx',4);

Rg = 8.314;
time = sheet4(:,1);				% mins
Texp = sheet4(:,2)+273.15;				% Celsius
Tout = sheet4(:,3)+273.15;				% Celsius
O2_out = sheet4(:,4);			% percentage vol.
CH4_out = sheet4(:,5);			% percentage vol.
C0 = [4 12 0 0 84]	;
C0 = [C0 C0]'
Cexp = [CH4_out O2_out -(CH4_out - C0(1)) 2*-(CH4_out - C0(1)) 100-(CH4_out+O2_out+3*-(CH4_out-C0(1)))];
par0 = randn(1,2);
			% percentage vol.
nu =[-1 -2 1 2 0]';
set(0,'DefaultAxesColorOrder',[1 0 0; .1 1 .1; 0 0 1; .1 0 .1; 0 1 .1])


plot(time,Cexp,'s')
ylim([0 15])
opt = optimset('Display','Iter');



for u=1:length(Texp)
T = Texp(u);

[par,fval] = fminsearch(@err,par0,opt,time,Cexp,Rg,T,a,hm);
parmatrix=[par0;par];


end

disp(parmatrix);
	function [S] = err(par,time,Cexp,Rg,T,a,hm)
		disp(par);
	A = par(1);
	Ea = par(2);
	k = A.*exp(-Ea./(Rg.*T));
	sol =solution(time,par)
	Ccalc = deval(sol,time)';
	S = norm(Ccalc - Cexp);
	t1=linspace(min(time),max(time),100);
	C1=deval(sol,t1);
	
	plot(t1,C1,time,Cexp,'s');
	ylim([0 15]);
	drawnow
	end
function [sol]=solution(time,par)

	C0 = [4 12 0 0 84];
	C0 = [C0 C0]'
	nu =[-1 -2 1 2 0]';
	Rg = 8.314;

	M = [1 0
     0 0];
options = odeset('Mass',M);	
sol = ode15s(@Bmi,time,C0,options,par,nu,Rg,T,a,hm);

end

function  [ADE]=Bmi(time,Cexp,par,nu,Rg,T,a,hm)
disp(par)
disp(nu)
disp(Rg)
Cexp1 = Cexp(1:5)
Cexp2 = Cexp(6:10)
	A = par(1);
	Ea = par(2);
	k = A.*exp(-Ea./(Rg.*T));
	kapp = k;
	%kapp = (1/mtc*ssa + 1/k)^-1;
R = kapp.*Cexp1(1)
r = nu * R


ADE = [    -a*hm * Cexp1
              hm * Cexp2+r  ]
end

end