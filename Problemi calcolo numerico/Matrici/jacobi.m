function[x,k]=jacobi(a,b,x0,toll,kmax)
k=0;
x=x0;
test=toll+1;
n=size(a,1);

while test>=toll&&k<kmax
    k=k+1;
    for i=1:n
        xn(i)=0;
        for j=1:i
            xn(i)=xn(i)+a(i,j)*xn(j);
        end
        for j=(i+1):n
            xn(i)=xn(i)+a(i,j)*xn(j);
        end
        xn(i)=(b(i)-xn(i))/a(i,i);
    end
    test=abs(xn-x);
    x=xn;
end

