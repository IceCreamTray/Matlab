function cpi=calspec(Tadim,Trapp,C)                  %calcolo dei calori specifici delle specie
Tdim=Tadim*Trapp-273.15;                             %Temperatura espressa in celsius
temp=[1 Tdim Tdim^2 Tdim^3];       
cpi=temp*C;                                          %calori specifici