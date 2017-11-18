n=50;
a=diag(4*ones(n,1))+diag(-2*ones(n-1,1),1)+diag(ones(n-1,1),-1);
sol=ones(n,1);
b=a*sol;
kmax=100;
[x,k]=jacobi(a,b,x0,toll,kmax);