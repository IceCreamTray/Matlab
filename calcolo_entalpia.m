function h=entalpia (Tadim,Trapp,C,Dhf)                                          %function che consente il calcolo delle entalpie di ogni specie
Tdim=Tadim*Trapp-273.15;                                                         %temperatura dimensionale in celsius
Tr=25;                                                                           
temp=[Tdim-Tr, 1/2*(Tdim^2-Tr^2), 1/3*(Tdim^3-Tr^3), 1/4*(Tdim^4-Tr^4)];         %calcolo dell'integrale
intcp=temp*C;                                                                    
h=Dhf+intcp;                                                                     %calcolo dell'entalpia                                                             