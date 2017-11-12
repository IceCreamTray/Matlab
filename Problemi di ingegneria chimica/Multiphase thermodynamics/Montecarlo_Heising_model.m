clear all, close all, clc
% Imposto il seed della matrice: genero sempre la stessa matrice casuale
s=rng;
A=rand(22,22);

rng(s);
% % Matrice su cui applico Montecarlo
MC=rand(22,22);

% Cambio i valori nella matrice in 1 e -1
for i=1:22
    for j=1:22
    if MC(i,j)<0.5
        MC(i,j)=-1;
    elseif MC(i,j)>=0.5
        MC(i,j)=1;
    end
    end
end

% Impongo che il contorno della matrice sia equivalente alla  riga opposta
% della matrice 20x20
MC(1,:)=MC(21,:);
MC(22,:)=MC(2,:);
MC(:,1)=MC(:,21);
MC(:,22)=MC(:,2);


format short
%disp('Matrice random per Montecarlo: '), disp(MC)
figure
pcolor(MC)

% Controllo che l'uguglianza dei bordi
%disp('Si verifica uguaglianza tra i vettori limite della matrice');
%disp('1: '), disp(MC(1,:)), disp(MC(21,:));
%disp('2: '), disp(MC(22,:)), disp(MC(2,:));
%disp('3: '), disp(MC(:,1)), disp(MC(:,21));
%disp('4: '), disp(MC(:,22)), disp(MC(:,2));


%% MONTECARLO
Kb = 1.38064852e-23;
Tr = 1;
T=600;
beta = 1/(T*Kb);


iterations = 100000;
interval = 100;

count=0;
index = 0;
vCount = zeros(1,interval);
vMeans =zeros(1,ceil(iterations/interval));
indexvec=zeros(1,ceil(iterations/interval));


for k=1:iterations
  
        % Generate a random couple ij to change in the matrix
        a=randi([1 21]); % x
        b=randi([1 21]); % y
        MC(a,b)=MC(a,b)*-1;
        %controls that the index are not out of boundary
        x_previous = a - 1;
        if x_previous < 1
            x_previous = 22;
        end

        x_next = a + 1;
        if x_next > 22
            x_next = 1;
        end

        y_previous = b - 1;
        if y_previous < 1
            y_previous = 22;
        end

        y_next = b + 1;
        if y_next > 22
            y_next = 1;
        end
        %Calculate the energy of the change as function of Tr
        DeltaEv=-2*1/2.269*Tr/T*Kb*(MC(a,b)*MC(a,y_next)+MC(a,b)*MC(a,y_previous)+MC(a,b)*MC(x_next,b)+MC(a,b)*MC(x_previous,b));
        
        
        if DeltaEv <= 0
            MC(a,b)=MC(a,b);
        else 
            if exp(DeltaEv*-beta)>= rand([0 1])
                MC(a,b)=MC(a,b);
            else
                MC(a,b)=MC(a,b)*-1;
                
                copy = 0;
                to_a = a;
                to_b = b;
                if a < 3
                    copy = 1;
                    to_a = a + 20;
                elseif a > 20
                    copy = 1;
                    to_a = a - 20;
                end
                if b < 3
                    copy = 1;
                    to_b = b + 20;
                elseif b > 20
                    copy = 1;
                    to_b = b - 20;
                end
                if copy == 1
                   MC(to_a,to_b) = MC(a,b); 
                end
            end
        end
    
   % Store the mean of DeltaEv for every [interval] iterations.
   
   count = count + 1;
   vCount(count) = DeltaEv;
   
   if count >= interval
       count = 0;
       
       index = index + 1;
       indexvec(index) = k;
       vMeans(index) = mean(vCount);
   end
 
end

% Store remaining stuff.
if count > 0
    index = index + 1
    indexvec(index) = iterations; % or does k remain?
    vMeans(index) = mean(vCount);
end

figure
pcolor(MC)
figure
plot(indexvec,vMeans)
disp(vMeans)

    
    
