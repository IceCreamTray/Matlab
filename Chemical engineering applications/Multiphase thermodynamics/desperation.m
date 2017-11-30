%% ====== Montecarlo simulation per modello bidimensionale di Ising ======
% Crea una matrice 20x20 con direzioni di spin randomizzate. Cambia un
% valore di spin per ogni iterazione e calcola la variazione di energia.
% Accetta o rifiuta il cambiamento in base al cambio di energia del
% sistema, secondo le condizioni di Montecarlo. Plotta l'energia finale del
% sistema e la varianza di energia.

%% ============================= MAIN SCRIPT =============================

clear all,clc, close all
addpath('../../util');			% Setta il percorso delle funzioni

%% Costanti operative
const_xdim = 22;									% Numero righe
const_ydim = 22;									% Numero colonne
const_iterations = 1e6;								% Numero di iterazioni
const_interval = 1e3;                               % Ampiezza intervallo

%% Variabili fisiche
% Temperatura ridotta
Tr = 0.5;
kb = physconst('Boltzmann');

% Costante di accoppiamento
J = kb / (2.269 * Tr);									% Adimensionale
beta = 1 / kb;											% Adimensionale 
betaj = beta * J;

%% Pre-allocamento
MC = zeros(22,22);
VEmean_buffer = zeros(1, const_interval);
VEmean = zeros(1, ceil(const_iterations / const_interval));
m_buffer = zeros(1, const_interval);
m_mean = zeros(1, ceil(const_iterations / const_interval));
var = zeros(1,ceil(const_iterations/const_interval));

%% Matrice
MC = random_bound_matrix(const_xdim, const_ydim);

%% Plot iniziali
figure
pcolor(MC);
title('Sistema iniziale');

%% Energia iniziale
E0_sys = calculate_system_energy(MC, const_xdim, const_ydim, 0.5); 
% Calcola l'energia iniziale del sistema chiamando la rispettiva funzione

%% Primary loop
% Setta i contatori
count = 0;
index = 0;

% Iterazioni
for i = 1:const_iterations
	
	%% CAMBIO
	% Cambia uno spin alle coordinate random a,b scelte nella matrice 20x20
	a = randi([2 (const_xdim - 1)]);
	b = randi([2 (const_ydim - 1)]);
	MC(a,b)= -MC(a,b);
	
	% Cambia le condizioni al contorno chiamando la rispettiva funzione
	MC = update_bound_matrix_bounds(MC, const_xdim, const_ydim);
	
	% Calcola l'energia dopo il cambio di spin
	Ec_sys = calculate_system_energy(MC, const_xdim, const_ydim, 0.5);
	
	% Calcola la variazione di energia
	DeltaE = Ec_sys - E0_sys;
	
	%% CONDIZIONI DI MONTECARLO
	
	if DeltaE > 0 		
		if exp(-betaj*DeltaE) < rand
			MC(a,b) = -MC(a,b);				
			MC = update_bound_matrix_bounds(MC, const_xdim, const_ydim);
			% se le condizioni non sono rispettate riporta lo spin al suo
			% valore iniziale e aggiorna il contorno
		end
	end
	
	% Calcola l'energia finale del sistema
	Efin_sys = calculate_system_energy(MC, const_xdim, const_ydim, 0.5);
	m = magnetization(MC,const_xdim,const_ydim);
	% L'energia finale � l'energia iniziale della successiva iterazione
	E0_sys = Efin_sys;
    

	%% RACCOLTA
	count = count + 1;								% Incrementa il contatore
	VEmean_buffer(count) = Efin_sys;				% Il buffer contiene tanti valori quanto ampio � l'intervallo desiderato
    m_buffer(count) = m;
	if count == const_interval
		count = 0;									% Quando il contatore arriva all'ampiezza desiderata
		index = index + 1;							% viene resettato e viene calcolata la media dei valori e
		VEmean(index) = mean(VEmean_buffer);		% viene riempito un vettore con la media di ogni intervallo
        m_mean(index) = mean(m_buffer);
	end
	
end

%% Calcola varianza e Cv
var = (VEmean - VEmean(end)).^2;
cv = (var/(2.269^2*Tr^2));

%% PLOT FINALI
% Matrice
figure
title('Sistema finale');
pcolor(MC)	

% Grafici
figure
plot(1:index,VEmean, '-r');
title('Energia del sistema');
ylabel('Energia media [-]');
xlabel('Intervalli');

figure
plot(1:index,m_mean,'-r');
title('Magnetizzazione');
ylabel('Magnetizzazione media [-]');
xlabel('Intervalli');

figure
plot(1:index,var,'-r');
title('Varianza');
xlabel('Intervalli');
ylim([-1e4 0.5*1e5]) 

figure
plot(1:index,cv,'-r');
title('Cv');
xlabel('Intervalli');
ylim([0 8e4]);
