%This program can be used for optical modeling of NbN devices with varying 
%geometries and varying wavelengths with wavelength-dependent materials.

%The code assumes back illumination with an incoherent source.

clear all
%%
%Geometric constants
d_ARC = 252;
d_SiNx = 252;
d_NbN = 4;
d_NbO = 0;
k_NbO = 1;
d_HSQ = 0;%264;
d_gold = 0;
fill_factor = 0.4;
n_air = 1; 

%%
%Initialize values for result matrices using number of wavelengths
max=476;
T_plot = zeros(max,1);
R_plot = zeros(max,1);
A_plot = zeros(max,1); 
A_NbN_plot = zeros(max,1); 
wavelength = zeros(max,1);

%%
%Read material values from csv file
%M = |wavelength|nNbN|kNbN|nSiO2|
M = csvread('NbN_data.csv');
gold_input = csvread('gold_input.csv');
Si_input = csvread('Si_input.csv'); %From Palik via refractiveindex.info

%%
for count=1:1:max
    
    %Load material constants
    lambda = M(count,1);
    k = 2*pi/lambda;
    k_air = n_air*k;
    
    n_ARC = M(count,4);
    k_ARC = n_ARC*k;
    
    n_substrate = Si_input(count,2) - 1i*Si_input(count,3);
    k_substrate = n_substrate*k;
    
    n_SiNx = n_ARC;
    k_SiNx = k_ARC;
    
    n_NbN_given = M(count,2) - 1i*M(count,3);
    epsilon_weighted_NbN = fill_factor*(n_NbN_given^2)+(1-fill_factor)*(n_air^2);
    n_NbN = sqrt(epsilon_weighted_NbN);
    k_NbN = n_NbN*k;
    
    % This value just an estimate
    n_HSQ = M(count,4)*0.9;
    k_HSQ = n_HSQ*k;
    
    n_gold =  gold_input(count,2) - 1i*gold_input(count,3);
    k_gold = n_gold*k;

    %Initialize values for loop
    T_total = 0;
    R_total = 0;
    A_total = 0;
    A_NbN_total = 0;
    inc = 1; %fraction of incident power
    
    %At the first interface (back of the sample)
    D_air_ARC = 0.5*[1+(k_ARC/k_air) 1-(k_ARC/k_air); 1-(k_ARC/k_air) 1+(k_ARC/k_air)];
    P_ARC = [exp(1i*k_ARC*d_ARC) 0; 0 exp(-1i*k_ARC*d_ARC)];
    D_ARC_substrate = 0.5*[1+(k_substrate/k_ARC) 1-(k_substrate/k_ARC); 1-(k_substrate/k_ARC) 1+(k_substrate/k_ARC)];
    M_ARC = D_air_ARC*P_ARC*D_ARC_substrate;
    
    R_ARC = abs(M_ARC(2,1)/M_ARC(1,1))^2;
    R_total = R_total+inc*R_ARC;
    T_ARC = (n_substrate/n_air)*abs(1/M_ARC(1,1))^2;
    inc = inc*T_ARC;
    
    %start for loop for muliple passes within the substrate
    for pass=1:1:10%number of internal reflections; generally three is adequate, but the calculation is fast
    
        %At the second interface (substrate/SiNx/NbN/cavity)
        D_substrate_SiNx = 0.5*[1+(k_SiNx/k_substrate) 1-(k_SiNx/k_substrate); 1-(k_SiNx/k_substrate) 1+(k_SiNx/k_substrate)];
        P_SiNx = [exp(1i*k_SiNx*d_SiNx) 0; 0 exp(-1i*k_SiNx*d_SiNx)];
        D_SiNx_NbN = 0.5*[1+(k_NbN/k_SiNx) 1-(k_NbN/k_SiNx); 1-(k_NbN/k_SiNx) 1+(k_NbN/k_SiNx)];
        P_NbN = [exp(1i*k_NbN*d_NbN) 0; 0 exp(-1i*k_NbN*d_NbN)];
        D_NbN_NbO = 0.5*[1+(k_NbO/k_NbN) 1-(k_NbO/k_NbN); 1-(k_NbO/k_NbN) 1+(k_NbO/k_NbN)];
        P_NbO = [exp(1i*k_NbO*d_NbO) 0; 0 exp(-1i*k_NbO*d_NbO)];
        D_NbO_HSQ = 0.5*[1+(k_HSQ/k_NbO) 1-(k_HSQ/k_NbO); 1-(k_HSQ/k_NbO) 1+(k_HSQ/k_NbO)];
        P_HSQ = [exp(1i*k_HSQ*d_HSQ) 0; 0 exp(-1i*k_HSQ*d_HSQ)];
        D_HSQ_gold = 0.5*[1+(k_gold/k_HSQ) 1-(k_gold/k_HSQ); 1-(k_gold/k_HSQ) 1+(k_gold/k_HSQ)];
        P_gold = [exp(1i*k_gold*d_gold) 0; 0 exp(-1i*k_gold*d_gold)];
        D_gold_air = 0.5*[1+(k_air/k_gold) 1-(k_air/k_gold); 1-(k_air/k_gold) 1+(k_air/k_gold)];
        M_second_interface = D_substrate_SiNx*P_SiNx*D_SiNx_NbN*P_NbN*D_NbN_NbO*P_NbO*D_NbO_HSQ*P_HSQ*D_HSQ_gold*P_gold*D_gold_air;

        R_second_interface = abs(M_second_interface(2,1)/M_second_interface(1,1))^2;
        T_second_interface = (n_air/n_substrate)*abs(1/M_second_interface(1,1))^2;
        T_total = T_total+inc*T_second_interface;
        A_total = A_total+inc*(1-T_second_interface-R_second_interface);
        
        %Absorption in NbN
        matrix = D_NbN_NbO*P_NbO*D_NbO_HSQ*P_HSQ*D_HSQ_gold*P_gold*D_gold_air;
        AB = matrix*[1/M_second_interface(1,1) 0]'; %gives the constants for the electric field in NbN
        imeps = (8.85418782*10^-12)*imag(n_NbN^2);
        c = 3*10^8; %c in m/s
        omega = 2*pi*c/(lambda*10^-9);
        d = d_NbN*10^-9; %d is in meters
        X = 0:d/1000:d;
        k_NbNm = k_NbN*10^9; %k_NbN in meters, not nm
        Y = abs(AB(1)*exp(-1i*k_NbNm*(X-d))+AB(2)*exp(1i*k_NbNm*(X-d))).^2;
        Z = inc*trapz(X,Y);
        irradiance = n_substrate*(3*10^8)*(8.85418782*10^-12)*0.5;%need to include n_substrate for agreement with previous simulations
        Q = 0.5*omega*imeps*Z/irradiance;
        A_NbN_total = A_NbN_total - Q;
        
        inc = inc*R_second_interface;
    
        %At the first interface
        D_substrate_ARC = 0.5*[1+(k_ARC/k_substrate) 1-(k_ARC/k_substrate); 1-(k_ARC/k_substrate) 1+(k_ARC/k_substrate)];
        P_ARC = [exp(1i*k_ARC*d_ARC) 0; 0 exp(-1i*k_ARC*d_ARC)];
        D_ARC_air = 0.5*[1+(k_air/k_ARC) 1-(k_air/k_ARC); 1-(k_air/k_ARC) 1+(k_air/k_ARC)];
        M_first_interface = D_substrate_ARC*P_ARC*D_ARC_air;
    
        R_first_interface = abs(M_first_interface(2,1)/M_first_interface(1,1))^2;
        T_first_interface = (n_air/n_substrate)*abs(1/M_first_interface(1,1))^2;
        R_total = R_total + inc*T_first_interface;
        inc = inc*R_first_interface;
    end
    
    %These values are in fractions, not percentages.
    T_plot(count)=T_total;
    R_plot(count) = R_total;
    A_plot(count) = A_total;
    A_NbN_plot(count) = A_NbN_total;
    wavelength(count)=lambda;
    R_plot+T_plot+A_plot; %check to make sure they add up to 100%
end


%%
plot(wavelength, A_NbN_plot*100, 'b-')
hold on

xlabel('wavelength (nm)')
ylabel('absorptance in NbN (%)')
