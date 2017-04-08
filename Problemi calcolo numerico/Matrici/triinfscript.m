%% TRIINFSCRIPT

%inizializzazione
L=[3.5 0 0; 4.04 26.4 0; 6.08 37.8 100.6];
b=[2.46; 11.55; 31.1];

%calcolo dimensione
n=size(L,1);

%verifica che sia triangolare superiore
if nnz(triu(L,1))==0
    [xv]=triinf(L,b(:));
else
    disp('la matrice non è triangolare inferiore');
    
end

disp('matrice:');
disp(L);
disp('vettore termine noto:');
disp(b);
disp('dimensione sistema');
disp(n);
disp('vettore soluzione');
disp(xv);