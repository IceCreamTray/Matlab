%% ========= Montecarlo simulation for bidimensional Ising model =========
% 

 close all, clear all, clc

% Imposto il seed della matrice: genero sempre la stessa matrice casuale
% s=rng;
% A=rand(22,22);
% 
% rng(s)
% Montecarlo matrix
MC = rand(22,22);

% Cambio i valori nella matrice in 1 e -1
for i = 1:22
    for j = 1:22
        if MC(i,j) < 0.5
            MC(i,j) = -1;
        elseif MC(i,j) >= 0.5
            MC(i,j) = 1;
        end
    end
end

% Boundary conditions
MC(1,:) = MC(21,:);
MC(22,:) = MC(2,:);
MC(:,1) = MC(:,21);
MC(:,22) = MC(:,2);

% Reference Temperatures
Tr = 0.5;      %[-]
Tr2 = 1;        %[-]
Tr3 = 1.5;      %[-]

% Adimensional Constants
J = 1 / 2.269 / Tr;    % It's J/Tkb
beta = 1;          % It's beta*Tkb


%%
iterations=10^6;
E0 = zeros(22,22);
E1=zeros(22,22);
pause = 500;
vectE = zeros(1,pause);
vectdeltaE = zeros(1,pause);
meanE = zeros(1,ceil(iterations/pause));
meandeltaE = zeros(1,ceil(iterations/pause));
count = 0;
index = 0;
%%
pcolor(MC)

% Total initial energy


for i=1:iterations

for u=2:21
    for k=2:21
        E0(u,k) = 0.5*(-J * (MC(u,k) * MC(u-1,k) + MC(u,k) * MC(u,k+1) + MC(u,k) * MC(u+1,k) + MC(u,k) * MC(u,k-1)));
    end
end
E0_sys = sum(sum(E0)');

a = randi([2 21]);
b = randi([2 21]);
MC(a,b) = -MC(a,b);

% conditions

MC(1,:) = MC(21,:);
MC(22,:) = MC(2,:);
MC(:,1) = MC(:,21);
MC(:,22) = MC(:,2);
 
    
 for u=2:21
     for k=2:21
        E1(u,k) = (-J * (MC(u,k) * MC(u-1,k) + MC(u,k) * MC(u,k+1) + MC(u,k) * MC(u+1,k) + MC(u,k) * MC(u,k-1)))*0.5;
     end
 end
 
 E1_sys = sum(sum(E1)');
 
 % DElta
 deltaE = E1_sys - E0_sys;
 
 % Conditions
 x = rand(0,1);
    if deltaE <= 0
        MC(a,b)=MC(a,b);
    else 
        if exp(-beta * deltaE) >= x
            MC(a,b)=MC(a,b);
        else         
            MC(a,b) = -MC(a,b);   % The spin change in neglected, so it is brought to its initial direction
MC(1,:) = MC(21,:);
MC(22,:) = MC(2,:);
MC(:,1) = MC(:,21);
MC(:,22) = MC(:,2);
 end
        end
     
    
 
   for u=2:21
     for k=2:21
         Ef(u,k) = (-J * (MC(u,k) * MC(u-1,k) + MC(u,k) * MC(u,k+1) + MC(u,k) * MC(u+1,k) + MC(u,k) * MC(u,k-1)))*0.5;
     end
   end  
 
   Ef_sys = sum(sum(Ef)');
   
   deltaE = Ef_sys - E0_sys;

   count = count + 1;
   vectE(count) = Ef_sys;
   vectdeltaE(count) = deltaE;
% Interactions Vectors
    if count>pause
        count=0;
        index = index + 1;
        meanE(index) = mean(vectE);
        meandeltaE(index) = mean(vectdeltaE);
    end
if count > 0
    index = index + 1;
    meanE(index) = mean(vectE);
    meandeltaE(index) = mean(vectdeltaE);
end
end

c=(1:index);
figure
plot(c,meanE), title('Energy of the system')
figure
plot(c,meandeltaE), title('Differences of Energy of the System')
figure
pcolor(MC), title('Final State')
