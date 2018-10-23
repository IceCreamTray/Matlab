%% Evaluation of experimental uncertainty with GUM and SUP methods.
% The method requires two measures per sample.
% The evaluation of the uncertainty interval is carried out with both SUP
% and GUM methods, based on the pseudocode retrievable in the article
% "Evaluating methods of calculating measurement uncertainty" by B.D.Hall,
% Measurement Standards Laboratory of New Zealand, Industrial Research Ltd,
% Lower Hutt, New Zealand.
% In case of low number of samples it is advised to rely on the SUP
% implementation, since the GUM has phisical meaning for a large number of samples.
% The methods are checked  to be reliable by counters at the end of the
% code.


clc, clear all, close all




%% Data
% Reads from Excel file

data =  xlsread("misure.xlsx");		% Excel file containing measurements.

l = data(end,1);

% Trims the necessary columns.
L = data(:, 1);
z1 = data(:, 2);					
z2 = data(:, 3);

mes = table(L, z1, z2, 'Variablenames', { 'Sample', 'Measure_1', 'Measure_2' });

disp('Measures table:');
disp(mes);




%% Counters inizialization
% The counters help evaluating the methods' reliability.

gamma_counter_gum = 0;
gamma_counter_sup = 0;



%% Variables

z = (z1+z2)./2;		% Calculates z.
u = std(z);			% Calculates the standard deviation.

% Calculates the means of both measurements vectors.
gamma1 = mean(z1);
gamma2 = mean(z2);




%% Sorting

zsorted = sort(z);				% Sorts z in ascending order.

Lup = l .* 0.975;				% Calculates upper index.
Lupround = round(Lup);			% Rounds upper index to the closest integer.

Llow = l .* 0.025;				% Calculates lower index.
Llowround = round(Llow);		% Rounds lower index to the closest integer.

if Llowround == 0				% Rounds the lower values to the closest integer different from zero.
	Llowround = 1;
end
	



%% Vectors

up = zeros(1,length(z))';
down = zeros(1,length(z))';




%% GUM METHOD
% The result is an interval in which thereis 95% of possibility to find the
% measured value. Outside said interval the mean is not reliable anymore.

for (k = 1:length(z))
	
	i = (z(k) - z1(k)) / z2(k);					% Calculates the imaginary part of the number for each sample.
	gamma = gamma1 + i * gamma2;				% Calculates the mean for the sample by combining real and imaginary part.

	% Calculates upper and lower values of the error interval by definition.
	down(k) = z(k) - 1.96 .* sqrt(u);
	up(k) = z(k) + 1.96 .* sqrt(u);

	if (gamma <= up(k) && gamma >= down(k))		% Checks the presence of the mean in said interval.
		gamma_counter_gum = gamma_counter_gum + 1;
	end
	
end




%% SUP METHOD

up_sup = zsorted(Lupround);			% Selects the upper bound of the interval as the z value related to the lowest upper index.
down_sup = zsorted(Llowround);		% Selects the lower bound of the interval as the z value related to the lowest bound index.
	

for (p = 1:length(z))
	
	i = (z(p) - z1(p)) / z2(p);					% Calculates the imaginary part of the number for each sample.
	gamma = gamma1 + i * gamma2;				% Calculates the mean for the sample by combining real and imaginary part.
	
	if (gamma <= up_sup && gamma >= down_sup)	% Checks the presence of the mean in said interval.
		gamma_counter_sup = gamma_counter_sup + 1;
	end
	
end
	



%% Checks reliability and displays results

disp('Reliability check: number of times the mean falls into the interval');
disp(' ');
rel = table(gamma_counter_gum,gamma_counter_sup,'Variablenames',{'GUM', 'SUP'});
disp(rel);


if (gamma_counter_gum < gamma_counter_sup)
	
	disp('WARNING: GUM is not reliable!');
	disp(' ');
	disp('Confidence interval with 95% certainty (SUP METHOD)');
	disp(' ');
	tab = table(down_sup',up_sup','Variablenames',{'lower_bound', 'upper_bound'});
	disp(tab);
	
elseif (gamma_counter_gum == gamma_counter_sup)
	
	disp('Both methods are reliable');
	disp(' ');
	disp('Confidence interval with 95% certainty (GUM METHOD)');
	disp(' ');
	T = table(mean(down),mean(up),'Variablenames',{'lower_bound', 'upper_bound'});
	disp(T);
	
	disp(' ');
	disp('Confidence interval with 95% certainty (SUP METHOD)');
	disp(' ');
	tab = table(down_sup',up_sup','Variablenames',{'lower_bound', 'upper_bound'});
	disp(tab);	
	
else
	
	disp('WARNING: SUP is not reliable!');
	disp(' ');
	disp('Confidence interval with 95% certainty (GUM METHOD)');
	disp(' ');
	T = table(mean(down),mean(up),'Variablenames',{'lower_bound', 'upper_bound'});
	disp(T);
	
end


s = (z1 + z2)./2;
disp('Values mean: '), disp(mean(s));
disp('St. deviation: '), disp(std((z1+z2)./2));









