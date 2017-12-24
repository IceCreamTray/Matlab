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
sheet1 = xlsread('Expdata.xlsx',1);
%sheet2 = xlsread('Expdata.xlsx',2);
%sheet3 = xlsread('Expdata.xlsx',3);

Rg = 8.314;
time = sheet1(:,1);				% mins
Tin = sheet1(:,2)+273.15;				% Celsius
Tout = sheet1(:,3)+273.15;				% Celsius
O2_out = sheet1(:,4).*(0.5/100);			% percentage vol.
CH4_out = sheet1(:,5).*(0.5/100);			% percentage vol.
C0 = [0.005 0.01 0 0 0.485]	;


Cexp = [CH4_out O2_out (CH4_out - C0(1)) 2*(CH4_out - C0(1)) (0.5-(CH4_out+O2_out+3*(CH4_out-C0(1))))];
k0 = 0.00001;
			% percentage vol.
nu =[-1 -2 1 2 0]';
set(0,'DefaultAxesColorOrder',[1 0 0; .1 1 .1; 0 0 1; .1 0 .1; 0 1 .1])

plot(time,Cexp,'s');
opt = optimset('Display','Iter');
[par,fval] = fminsearch(@err,k0,opt,time,Cexp);
disp(par);



function [S] = err(par,time,Cexp)
sol =solution(time,par)
Ccalc = deval(sol,time)';
S = norm(Ccalc - Cexp);
t1=linspace(min(time),max(time),100);
C1=deval(sol,t1);

plot(t1,C1,time,Cexp,'s');
ylim([-0.02 0.02]);
drawnow
end
function [sol]=solution(time,par)
	C0 = [0.005 0.01 0 0 0.485]	;
	nu =[-1 -2 1 2 0]';

	
sol = ode15s(@Bmi,time,C0,[],par,nu);
end
function  [PFR]=Bmi(par,Rg,Cexp,nu)


R = par*Cexp;
r = nu * R;
PFR = r;
end

end




