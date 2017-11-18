
%definizione matrice e termine noto
U=[1.75 2.02 3.04; 0 13.2 18.9; 0 0 50.3];
b=[4.92 23.1 62.2];

%calcolo dimensione matrice
dim=size(U,1);

%verifica che sia triangolare superiore
A=tril(U,-1);
if nnz(A)~=0
    disp('la matrice non è triangolare superiore!');
else
    xv =trisup(U,b);
    
    disp(U);
    disp(b);
    disp(dim);
    disp(xv);
end
