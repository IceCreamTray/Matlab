function CSTR_2_3_a
    %	Risoluzione dei bilanci di materia per un reattore CSTR
    %   NR = 2, NC = 4; R elementari
    %           2A  <-> B           reversibile!
    %           B  -> C+D
    %	ordine specie: (A B C D)

    % pulisce lo schermo
    close all, clear all, clc;
    options = optimset('Display','off');

    % J/molK
    R = 8.3145;

    vettconv = [];
    for th = 10:10:300
        for T = 5:1.89655:59.99995
            T = T + 273.15;

            % k1dir L/mol/min k2 1/min k3 1/min
            k = [0.3425*exp(-18800/(R*T)) 8.269e6*exp(-47830/(R*T))...
                18.056e9*exp(-79000/(R*T))];

            % matrice stechiometrica Ac DAA MO w
            nu=[-2  0
                1 -1
                0  1
                0  1];

            % concentrazioni in ingresso, Ac puro mol/L
            CIN = [13.5 0 0 0]';

            % stima della possibile soluzione
            C0  = [10 1 1 1]';
            C = fsolve(@BMi, C0, options, k, nu, CIN, th);

            % conversione di Ac
            X = 1 - C(1) / CIN(1);

            vettconv = [vettconv X];
        end
    end

    convfin=[];
    for i=1:1:30
        for j=1:1:30
            convfin(i,j)=vettconv((i-1)*30+j);
        end
    end


    [temp,tau]=meshgrid(5:1.89655:59.99995,10:10:300);
    surf(temp,tau,convfin);
    title('Conversion vs Temperature and Residence Time')
    ylabel('tau[s]');
    xlabel('T[°C]');
    zlabel('X[-]');
end

% ----------------------------------------------------------

function err = BMi(C,k,nu,CIN,th)
    % bilanci materiali sulle specie
    %	A  <-> 2B           reversibile!
    %	B  -> C
    % velocita' delle reazioni 1 e 2
    R  = [k(1)*C(1)^abs(nu(1,1)) - k(2)*C(2)^nu(2,1)
          k(3)*C(2)^abs(nu(2,2)) ];
    r  = nu*R;           % velocita' di produzione delle singole specie
    DC = C-CIN;
    err = DC - r*th ;     % 0 = BMi
end
