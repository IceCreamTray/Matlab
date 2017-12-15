clear all, clc, close all
addpath('../../util');								% Functions folder path

%% Operating variables
const_xdim = 22;									% Rows
const_ydim = 22;									% Columns
const_iterations = 1e5;								% iterations needed
const_interval = 1e2;                               % Sampling interval
equilibrium_sampling = 700;                         % Indicates the interval at which the equilibrium is reached

%% Physical variables
% Reduced temperature
Tr = 0.10

% Boltzmann constant
kb = physconst('Boltzmann');

%% Pre-allocation
MC = zeros(const_xdim,const_ydim);
VEmean_buffer = zeros(1, const_interval);
VEmean = zeros(1, ceil(const_iterations / const_interval));
m_buffer = zeros(1, const_interval);
m_mean = zeros(1, ceil(const_iterations / const_interval));
Trvec=zeros(1,10)
cv_vec=zeros(1,10);

for k=1:100
	Tr = Tr + 0.05;
	Trvec(k) = Tr;
	disp(k)
% Coupling constant
J = kb / (2.269 * Tr);									% Adimensional
beta = 1 / kb;											 
betaj = beta * J;										% Adimensional



%% Matrix
MC = random_bound_matrix(const_xdim, const_ydim);


%% Initial Energy
E0_sys = calculate_system_energy(MC, const_xdim, const_ydim, 0.5); 
% Calculates the initial energy by calling the proper function.

%% Initial Magnetization
m0 = magnetization(MC,const_xdim,const_ydim);
% Calculates the initial magnetization by calling the proper function.

%% Primary loop
% Counters setup
count = 0;
index = 0;

% Iterations
for i = 1:const_iterations
	
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
		if exp(-betaj*DeltaE) < rand
			MC(a,b) = -MC(a,b);				
			MC = update_bound_matrix_bounds(MC, const_xdim, const_ydim);
			% If the conditions are met it brings the spin back to its
			% initial direction and the boundaries are updated again. 
			% The variation of energy is too high.
		end
	end
	
	% Calculates the final energy of the system and the final
	% magnetization.
	Efin_sys = calculate_system_energy(MC, const_xdim, const_ydim, 0.5);
	m = magnetization(MC,const_xdim,const_ydim);
	
	% The final energy value is the initial value of the next iteration.
	E0_sys = Efin_sys;
    

	%% SAMPLING
	count = count + 1;								% Counter rising.
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
VEmeanstable = VEmean(equilibrium_sampling:end);						% By graphical analysis, it is observed that 
average = mean(VEmeanstable);						% the equilibrium is reached after 300 sampling intervals.
variance = var(VEmeanstable);
cv = (variance/(2.269^2*Tr^2));
cv_vec(k) = cv;


end
figure
scatter(Trvec,cv_vec);
title('Heat capacity at constant volume');
ylabel('Cv [-]');
xlabel('Tr [-]');
xlim([0 5]);