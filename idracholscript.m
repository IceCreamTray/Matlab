%% IDRACHOLSCRIPT
% definizione matrici
A=[-0.370 0.050 0.050 0.070; 0.050 -0.116 0 0.050; 0.050 0 -0.116 0.050; 0.070 0.050 0.050 -0.202];
b=[-2; 0; 0; 0];
n=size(A,1);

%trasformazione
%uso gauss,moltiplico tutto per -1
A1=A*-1;
b1=b*-1;

%cholesky
[L,flag]=chol(A1,'lower');

%calcolo L*Lt
M=L*L';
%%disp(M);
%%disp(A1);


%controllo
if flag==0
    %risolvo Ly=b
    [y]=triinf(L,b1);
    %risolvo L'x=y
    [xv]=trisup(L',y);
else
    disp('la matrice non è simmetrica e definita positiva');
    
end

disp('matrice')
disp(A)
disp('matrice equivalente')
disp(A1)
disp('vettore b')
disp('vettore equivalente')
disp(b1)
disp('dimensione del sistema')
disp(n)
disp('matrice L')
disp(L)
disp('soluzione ottenuta con cholensky')
disp(xv)
disp('soluzione ottenuta con matlab')
disp(A\b);
