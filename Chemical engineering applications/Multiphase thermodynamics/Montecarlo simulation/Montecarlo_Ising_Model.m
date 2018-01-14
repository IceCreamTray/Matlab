%% ====== Montecarlo simulation for bidimensional Ising Model ======
% Creates a 20x20 matrix with randomized spin directions. Changes a spin
% for each itaration and calculates the variation of energy in the system.
% Accepts or excludes the change of spin considering the Montecarlo
% conditions on the energy variation. Plots the energy and magnetization 
% profile, and the initial and final system configuration. 
% Gives out the values of variance and heat capacity at constant volume 
% when the equilibrium is reached.

%% ============================= MAIN SCRIPT =============================

%clear all, clc, close all
addpath('../../util');								% Functions folder path

%% Operating variables
const_xdim = 22;									% Rows
const_ydim = 22;									% Columns
const_iterations = 1e6;								% iterations needed
const_interval = 1e3;                               % Sampling interval
equilibrium_sampling = 700;                         % Indicates the interval at which the equilibrium is reached

%% Physical variables
% Reduced temperature
Tr = 0.5; 
%Tr = 1; 
%Tr = 1.5;

% Boltzmann constant
kb = physconst('Boltzmann');

% Coupling constant
J = kb / (2.269 * Tr);								% Nondimensional
beta = 1 / kb;											 
betaj = beta * J;									% Nondimensional

%% Pre-allocation
MC = zeros(const_xdim,const_ydim);
VEmean_buffer = zeros(1, const_interval);
VEmean = zeros(1, ceil(const_iterations / const_interval));
m_buffer = zeros(1, const_interval);
m_mean = zeros(1, ceil(const_iterations / const_interval));

%% Matrix
% Creates the matrix by calling the proper function.
MC = random_bound_matrix(const_xdim, const_ydim);

%% Initial plots
figure(1)
subplot(1, 2, 1);
pcolor(MC);
title('Initial system');

%% Initial Energy
% Calculates the initial energy by calling the proper function.
E0_sys = calculate_system_energy(MC, const_xdim, const_ydim, 0.5); 
disp('Initial Energy of the system: ');
disp(E0_sys);

%% Initial Magnetization
% Calculates the initial magnetization by calling the proper function.
m0 = magnetization(MC, const_xdim, const_ydim);
disp('Initial magnetization of the system: ');
disp(m0);

%% Primary loop
% Counters setup
count = 0;
index = 0;

% Iterations
for i = 1 : const_iterations
	
	%% Spin change
	% Changes a spin value at random coordinates picked from the 20x20
	% matrix.
	a = randi([2 (const_xdim - 1)]);
	b = randi([2 (const_ydim - 1)]);
	MC(a,b)= -MC(a,b);
	
	% Updates the boundaries by calling the proper function.
	MC = update_bound_matrix_bounds(MC, const_xdim, const_ydim);
	
	% Calculates the system energy after the spin change.
	Ec_sys = calculate_system_energy(MC, const_xdim, const_ydim, 0.5);
	
	% Calculates the energy variation.
	DeltaE = Ec_sys - E0_sys;
	
	%% MONTECARLO CONDITIONS
	
	if DeltaE > 0 		
		if exp(-betaj * DeltaE) < rand
			% If the conditions are met, the variation of energy is too
			% high and the spin is brought back to its initial direction.
			MC(a,b) = -MC(a,b);	
			
			% Updates the boundaries by calling the proper function.
			MC = update_bound_matrix_bounds(MC, const_xdim, const_ydim);
		end
	end
	
	% Calculates the final energy of the system and the final
	% magnetization.
	Efin_sys = calculate_system_energy(MC, const_xdim, const_ydim, 0.5);
	m = magnetization(MC, const_xdim, const_ydim);
	
	% The final energy value is the initial value of the next iteration.
	E0_sys = Efin_sys;
    

	%% SAMPLING
	count = count + 1;								
	VEmean_buffer(count) = Efin_sys;				% Buffers contain as many values as the size of the wanted interval.
    m_buffer(count) = m;							% They store both magnetization and energy values.
	
	if count == const_interval
		count = 0;									% When the counter reaches the sampling interval size
		index = index + 1;							% it is reset and the mean of the values in the buffer is calculated.
		VEmean(index) = mean(VEmean_buffer);		% A vector is filled with the mean values.
        m_mean(index) = mean(m_buffer);
	end
	
end

%% Variance and Cv calculation
VEmeanstable = VEmean(equilibrium_sampling : end);	% By graphical analysis, it is observed that 
magnetstable = m_mean(equilibrium_sampling : end);	% the equilibrium is reached after 300 sampling intervals.
average = mean(VEmeanstable);						% but to be sure to avoid fluctuations, samples are taken
variance = var(VEmeanstable);						% starting from the 700th interval.
cv = (variance / (2.269 ^ 2 * Tr ^ 2));
magnet = mean(magnetstable);

% Final Display
disp('Average system energy at equilibrium: ');
disp(average);
disp('Variance of the system at equilibrium: ');
disp(variance);
disp('Magnetization of the system at equilibrium: ');
disp(magnet);
disp('Heat capacity at constant volume at equilibrium: ');
disp(cv);

%% Final plots
whitebg([1 1 1]);

% Matrix
figure(1)
subplot(1, 2, 2);
pcolor(MC)
title('Final system');

% Plots
figure(2)
plot(1 : index, VEmean, 'LineWidth', 1.5), hold on
title('Adimensional system energy profile');
ylabel('Average energy on intervals [E/J]');
xlabel('Interval');
legend('Tr = 0.5', 'Tr = 1', 'Tr = 1.5');

figure(3)
loglog(1 : index, VEmean, 'LineWidth', 1.5), hold on
title('Adimensional system energy profile');
ylabel('Average energy on intervals [E/J]');
xlabel('Interval');
legend('Tr = 0.5', 'Tr = 1', 'Tr = 1.5');

figure(4)
plot(1 : index, m_mean, 'LineWidth', 1.5), hold on
title('Adimensional magnetization profile');
ylabel('Average magnetization on intervals [M/\mu]');
xlabel('Interval');
legend('Tr = 0.5', 'Tr = 1', 'Tr = 1.5');

figure(5)
semilogx(1 : index, m_mean, 'LineWidth', 1.5), hold on
title('Adimensional magnetization profile');
ylabel('Average magnetization on intervals [M/\mu]');
xlabel('Interval');
legend('Tr = 0.5', 'Tr = 1', 'Tr = 1.5');
