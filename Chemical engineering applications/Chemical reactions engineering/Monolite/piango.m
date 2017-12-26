clc,clear all
close all

%% parameters
% dati da gas
D   = 1e-5;   % m^2/s
L = 8.5e-3;      % lunghezza reattore, m
Sh = 3;
d = 8e-3;       % m, diametro canale
a = 4/d;        % 1/m  area specifica
R = 8.314;

% Read the experimental data from xls
%sheet1 = xlsread('Expdata.xlsx',1);
%sheet2 = xlsread('Expdata.xlsx',2);
%sheet3 = xlsread('Expdata.xlsx',3);
sheet4 = xlsread('Expdata.xlsx',4);


time = sheet4(:,1);				% mins
Texp = sheet4(:,2)+273.15;				% Celsius
Tout = sheet4(:,3)+273.15;				% Celsius
O2_out = sheet4(:,4);			% percentage vol.
CH4_out = sheet4(:,5);			% percentage vol.
C0 = [4 12 0 0 84]	;
Cexp = [CH4_out O2_out -(CH4_out - C0(1)) 2*-(CH4_out - C0(1)) 100-(CH4_out+O2_out+3*-(CH4_out-C0(1)))];
par0 = [1e10 1e5];

nu =[-1 -2 1 2 0]';
set(0,'DefaultAxesColorOrder',[1 0 0; .1 1 .1; 0 0 1; .1 0 .1; 0 1 .1])

plot(time,Cexp,'s')
ylim([0 15])
opt = optimset('Display','Iter');


T = Texp(u);
Creag = [C0(1) C0(1)]';
par=par0;
[t,Y] = cose(D,L,Sh,d,a,T,Creag,R,par,time)
disp(Y)
plot(t,Y)
