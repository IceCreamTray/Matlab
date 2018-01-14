clear all,close all

%% Read data
data = xlsread('matlab plots.xlsx',1);
data1 = xlsread('matlab plots.xlsx',2);
data2 = xlsread('matlab plots.xlsx',3);
data3 = xlsread('matlab plots.xlsx',4);
data4 = xlsread('matlab plots.xlsx',5);
data5 = xlsread('matlab plots.xlsx',6);

%% Methanol acetone at 328K
figure
x1_calc = data(:,4);
x1_exp = data(:,2);
pcalc = data(:,6);
pexp = data(:,1);
y1_calc = data(:,5);
y1_exp = data(:,3);
plot(x1_calc,pcalc,'LineWidth',1.5), hold on
plot(y1_calc,pcalc,'LineWidth',1.5), hold on
scatter(x1_exp,pexp),hold on
scatter(y1_exp,pexp);
title('Isothermal VLE - Pxy for an acetone(1)-methanol(2) mixture at 328K');
ylabel('P [kPa]');
xlabel('X_1 [mol/mol]')
legend('Calculated bubble curve','Calculated dew curve','x_1exp[-]','y_1exp[-]','location','best');

%% Methanol acetone at 373K
figure
x1_calc = data1(:,3);
x1_exp = data1(:,2);
pcalc = data1(:,6);
pexp = data1(:,1);
y1_calc = data1(:,4);
y1_exp = data1(:,5);
plot(x1_calc,pcalc,'LineWidth',1.5), hold on
plot(y1_calc,pcalc,'LineWidth',1.5), hold on
scatter(x1_exp,pexp),hold on
scatter(y1_exp,pexp);
title('Isothermal VLE - Pxy for an acetone(1)-methanol(2) mixture at 373K');
ylabel('P [kPa]');
xlabel('X_1 [mol/mol]')
legend('Calculated bubble curve','Calculated dew curve','x_1exp[-]','y_1exp[-]','location','south');

%% Water metanol at 323K
figure
x1_calc = data2(:,5);
x1_exp = data2(:,2);
pcalc = data2(:,6);
pexp = data2(:,1);
y1_calc = data2(:,4);
y1_exp = data2(:,3);
plot(x1_calc,pcalc,'LineWidth',1.5), hold on
plot(y1_calc,pcalc,'LineWidth',1.5), hold on
scatter(x1_exp,pexp),hold on
scatter(y1_exp,pexp);
title('Isothermal VLE - Pxy for a methanol(1)-water(2) mixture at 323K');
ylabel('P [kPa]');
xlabel('X_1 [mol/mol]')
legend('Calculated bubble curve','Calculated dew curve','x_1exp[-]','y_1exp[-]','location','best');

%% Water metanol at 353K
figure
x1_calc = data3(:,4);
x1_exp = data3(:,2);
pcalc = data3(:,6);
pexp = data3(:,1);
y1_calc = data3(:,5);
y1_exp = data3(:,3);
plot(x1_calc,pcalc,'LineWidth',1.5), hold on
plot(y1_calc,pcalc,'LineWidth',1.5), hold on
scatter(x1_exp,pexp),hold on
scatter(y1_exp,pexp);
title('Isothermal VLE - Pxy for a methanol(1)-water(2) mixture at 353K');
ylabel('P [kPa]');
xlabel('X_1 [mol/mol]')
legend('Calculated bubble curve','Calculated dew curve','x_1exp[-]','y_1exp[-]','location','best');

%% Water acetone at 318k
figure
x1_calc = data4(:,4);
x1_exp = data4(:,2);
pcalc = data4(:,6);
pexp = data4(:,1);
y1_calc = data4(:,5);
y1_exp = data4(:,3);
plot(x1_calc,pcalc,'LineWidth',1.5), hold on
plot(y1_calc,pcalc,'LineWidth',1.5), hold on
scatter(x1_exp,pexp),hold on
scatter(y1_exp,pexp);
title('Isothermal VLE - Pxy for an acetone(1)-water(2) mixture at 318K');
ylabel('P [kPa]');
xlabel('X_1 [mol/mol]')
legend('Calculated bubble curve','Calculated dew curve','x_1exp[-]','y_1exp[-]','location','best');

%% Water acetone at 373.15k
figure
x1_calc = data5(:,4);
x1_exp = data5(:,2);
pcalc = data5(:,6);
pexp = data5(:,1);
y1_calc = data5(:,5);
y1_exp = data5(:,3);
plot(x1_calc,pcalc,'LineWidth',1.5), hold on
plot(y1_calc,pcalc,'LineWidth',1.5), hold on
scatter(x1_exp,pexp),hold on
scatter(y1_exp,pexp);
title('Isothermal VLE - Pxy for an acetone(1)-water(2) mixture at 373K');
ylabel('P [kPa]');
xlabel('X_1 [mol/mol]')
legend('Calculated bubble curve','Calculated dew curve','x_1exp[-]','y_1exp[-]','location','best');


