function [xv] = triinf (L,b)

n=size(L,1);

if prod(diag(L))==0
    disp('il sistema è singolare!');
    
else

for i=1:n
    
    xv(i)=b;
    
    for j=1:(i-1)
        
        xv(i)=xv(i)-L(i,j)*xv(j);
        
    end
    
    xv(i)=xv(i)/L(i,i);
    
end
    
end