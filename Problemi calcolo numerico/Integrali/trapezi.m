function[int,h]= trapezi(f,a,b,m)
%TRAPEZI: metodo dei trapezi composto



%inizializzazione
h=(b-a)/m;                  %passo di integrazione
int= feval(f,a)+feval(f,b); %approssimazione dell'integrale come 
                            %somma della valutazione di f negli estremi dell'intervallo
%ciclo                            
for i=1:(m-1)
    x=a+i*h;                %formula dei trapezi composta
    int=int+2*feval(f,x);
end

%uscita
int=h*(int/2);              %integrale approssimato

end
