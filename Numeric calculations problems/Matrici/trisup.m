function [xv] = trisup (U,b)

n=size(U,1);
xv=zeros(n,1);

if prod(diag(U))==0
    disp('il sistema è singolare!');
    
else


for i=n:-1:1
    xv(i)=b(i);
    
    for j=(i+1):n
        
        xv(i)=xv(i)-U(i,j).*xv(i);
        
    end
    
    xv(i)=xv(i)./U(i,i);
    
end
    
end
    
        