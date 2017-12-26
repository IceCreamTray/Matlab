%% ================= CATALYTIC HONEYCOMB REACTOR =========================
% A small sample of a commercial catalyst for emission control in automotive 
% applications has been tested in the laboratory to measure its ability to 
% facilitate the oxidation of hydrocarbon. The most stable hydrocarbon (CH4) 
% has been selected as a benchmark of activity.
% The activity is measured in a flow reactor, at different temperature values.
% In addition to the qualitative indication of the minimum temperature required to
% activate some CH4 conversion, a tentative kinetic model is required, describing the surface 
% reaction, for scaling-up and muffler design.
% The catalyst sample used has a honeycomb structure, made of 23 square channels,
% 1x1x8.5 mm in size, fitting inside a quartz pipe, 8mm ID. 
% The temperature of the catalyst is monitored at the outlet while the inlet 
% temperature is used to control the oven. The inlet temperature is scanned from ambient to 800°C (850 in 2 cases),
% at a constant heating rate of 3°/min. The inlet volumetric flow rate is always 500 mLSTD/min.
%% 
function kool
%% ===================== MAIN =========================
clc,clear all,close all
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
ssa = 4/1; %1/mm
mtc = 6*10^-3*3*8;

Cexp = [CH4_out O2_out -(CH4_out - C0(1)) 2*-(CH4_out - C0(1)) 100-(CH4_out+O2_out+3*-(CH4_out-C0(1)))];
par0 = randn(1,2);
			% percentage vol.
nu =[-1 -2 1 2 0]';
set(0,'DefaultAxesColorOrder',[1 0 0; .1 1 .1; 0 0 1; .1 0 .1; 0 1 .1])


plot(time,Cexp,'s')
ylim([0 15])
opt = optimset('Display','Iter');

for u=1:length(Texp);
T = Texp(u);

[par,fval] = fminsearch(@err,par0,opt,time,Cexp,Rg,T);
parmatrix=[par0;par];


end
disp(parmatrix);
	function [S] = err(par,time,Cexp,Rg,T,mtc,ssa)
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
	nu =[-1 -2 1 2 0]';

	
sol = ode15s(@Bmi,time,C0,[],par,nu,Rg,T,mtc,ssa);

end
function  [PFR]=Bmi(time,Cexp,par,nu,Rg,T,mtc,ssa)
disp(par)
disp(nu)
disp(Rg)
	A = par(1);
	Ea = par(2);
	k = A.*exp(-Ea./(Rg.*T));
	kapp = k;
	%kapp = (1/mtc*ssa + 1/k)^-1;
R = kapp.*Cexp(1)
r = nu * R
PFR = r;
end

end




