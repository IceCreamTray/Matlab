%% Evaluation of experimental uncertainty
% The method requires two measures per sample.
% The evaluation of the uncertainty interval is carried out with both SUP
% and GUM methods, based on the pseudocode retrievable in the article
% "Evaluating methods of calculating measurement uncertainty" by B.D.Hall,
% Measurement Standards Laboratory of New Zealand, Industrial Research Ltd,
% Lower Hutt, New Zealand.
% In case of low number of samples it is advised to rely on the SUP
% implementation, since the GUM has phisical meaning for a large number of samples.

clc,  clear all, close all
%% Data
% Reads from Excel file
data =  xlsread("misure.xlsx");			% pre - filled excel file containing measurements.
disp('Measures table:');
z1 = data(:,2);							% Trims the necessary columns.						
z2 = data(:,3);
L = data(:,1);
mes = table(L,z1,z2,'Variablenames',{'Sample', 'Measure_1', 'Measure_2'});
disp(mes);

%% Counters inizialization
counter = 0;							% The counters help evaluating the methods efficiency.
coun = 0;

%% Variables 
z = sqrt(z1.^2 + z2.^2);				% Calculates the square error.
gamma1 = mean(z1);						% Calculates the means of both measurements vectors.
gamma2 = mean(z2);
u = var(z);								% Calculates the variance.

%% Sorting
zsorted = sort(z);						% Sorts z in ascending order.
Llow = L .* 0.025;						% Calculates lower indexes.
Llowround = zeros(1, length(Llow))';	% Helps storing the lower indexes.
Lup = L.*0.975;							% Calculates upper indexes.
Lupround = round(Lup);					% Rounds upper indexes to the closest integer.
for u = 1 : length(Llow)				% Rounds the lower values to the closest integer different from zero.
	
	Llowround(u) = round(Llow(u));
	
	if Llowround(u) == 0
		Llowround(u) = 1;
	end
	
end

%% Vectors
up = zeros(1,length(z))';
down = zeros(1,length(z))';

%% SUP METHOD
% The result is an interval in which thereis 95% of possibility to find the
% measured value. Outside said interval the mean is not reliable anymore.

for k = 1 : length(z1)
	
	i = (z(k) - z1(k)) / z2(k);				% Calculates the imaginary part of the number for each sample.
	gamma = gamma1 + i*gamma2;				% Calculates the mean for the sample by combining real and imaginary part.

	down(k) = z(k) - 1.96 .* sqrt(u);		% Calculates upper and bound values of the error interval by definition.
	up(k) = z(k) + 1.96 .* sqrt(u);

	if (gamma <= up(k) && gamma >= down(k))	% Checks the presence of the mean in said interval.
		counter = counter + 1;
	end
	
end

%% GUM METHOD
for g = 1:length(z)								
	
	up_gum(g) = zsorted(Lupround(g));				% Selects the upper bound of the interval as the z value related to the lowest upper index.
	down_gum(g) = zsorted(Llowround(g));			% Selects the lower bound of the interval as the z value related to the lowest bound index.
	
	if (gamma <= up_gum(g) && gamma >= down_gum(g))	% Checks the presence of the mean in said interval.
		coun = coun + 1;
	end
	
end

%% Checks reliability nad displays results
rel = table(counter,coun,'Variablenames',{'SUP', 'GUM'});
if counter < coun
	disp('WARNING: SUP is not reliable!');
	disp(' ');
	disp('Confidence interval with 95% certainty (GUM METHOD)');
	disp(' ');
	tab = table(down_gum',up_gum','Variablenames',{'lower_bound', 'upper_bound'});
	disp(tab);
	elseif counter == coun
		disp('Both methods are reliable');
		disp(' ');
		disp('Confidence interval with 95% certainty (SUP METHOD)');
		disp(' ');
		T = table(down,up,'Variablenames',{'lower_bound', 'upper_bound'});
		disp(T);
		disp(' ');
		disp('Confidence interval with 95% certainty (GUM METHOD)');
		disp(' ');
		tab = table(down_gum',up_gum','Variablenames',{'lower_bound', 'upper_bound'});
		disp(tab);
else
	disp('WARNING: GUM is not reliable!');
	disp(' ');
	disp('Confidence interval with 95% certainty (SUP METHOD)');
	disp(' ');
	T = table(down,up,'Variablenames',{'lower_bound', 'upper_bound'});
	disp(T);
end









