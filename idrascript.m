%% IDRASCRIPT
% inizializzazione 
A=[-0.370 0.050 0.050 0.070; 0.050 -0.116 0 0.050; 0.050 0 -0.116 0.050; 0.070 0.050 0.050 -0.202];
b=[-2;0;0;0];

[L,U,P]=lu(A);

%risolvo il sistema Ap=b,PxA=LxU
%LUx=Pb, Ly=Pb, Ux=y;

%Ly=Pb
[y]=triinf(L,P*b);

%Ux=y
[xv]=trisup(U,y);


disp('matrice A');
disp(A);
disp('vettore b');
disp(b);
disp('dimensione sistema');
disp(n);
disp(' vettore soluzione con metodi triangolari');
disp(xv);

disp('soluzione comando matlab')
disp(A\b);

