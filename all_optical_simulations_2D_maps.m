%This program can be used for optical modeling of NbN on a variety of
%devices. To exclude layers (e.g., the ARC), set the thickness of the layer
%in question to zero (e.g., d_ARC = 0).

%The code assumes back illumination with an incoherent source.

%The results are 2D plots of values.

clear all
%%
%Values and constants
lambda0 = 1550; %all distances are in nanometers
k = 2*pi/lambda0;

n_substrate = 3.47772; %MgO: 1.71, sapphire: 1.75, Si: 3.47772
k_substrate = n_substrate*k;

n_air = 1;
k_air = n_air*k; 

n_ARC = sqrt(n_air*n_substrate);%If there is no ARC but a layer of SiNx or SiO2 instead, input the values for the SiNx layer for n_ARC and d_ARC.
k_ARC = n_ARC*k;
d_ARC =lambda0/(4*n_ARC);

n_HSQ = 1.38;
k_HSQ = n_HSQ*k;
d_HSQ = 0; %optical cavity thickness

n_NbN_given = 5.23-5.82*1i; %Vikas' value
n_NbO_given = 2.28;
fill_factor = 0.5; %value between 0 and 1, where 1 is a continuous film

d_NbN =4;
d_NbO = 0;

n_SiNx=  2.217; %Silicon oxide at 1550 nm: 1.53626; SiNx: 2.217
k_SiNx = n_SiNx*k;
d_SiNx = lambda0/(4*n_SiNx); %red SiNx: 391, green SiNx: 300

n_gold = 0.55-11.5*1i;
k_gold = n_gold*k;
d_gold = 0;%d_gold is typically 120 nm for cavities on devices

%%
%Initialize values for result matrices
max_thickness=50;
max_fill=50;
T_plot = zeros(max_thickness+1,max_fill+1);
R_plot = zeros(max_thickness+1,max_fill+1);
A_plot = zeros(max_thickness+1,max_fill+1); %total absorption in structure, including absorption in gold if there is a gold layer
A_NbN_plot = zeros(max_thickness+1,max_fill+1);%absorption only in NbN layer
thickness = zeros(max_thickness+1, 1);
fill = zeros(max_fill+1, 1);

%%
for count_fill=0:1:max_fill
    fill_factor=count_fill/max_fill;
    %This calculation depends on the material between the nanowires, which
    %could be air or HSQ. Fill in n_air or n_HSQ in second term on rhs.
    epsilon_weighted = fill_factor*(n_NbN_given^2)+(1-fill_factor)*(n_air^2);
    n_NbN = sqrt(epsilon_weighted);
    epsilon_weighted = fill_factor*(n_NbO_given^2)+(1-fill_factor)*(n_air^2);
    n_NbO = sqrt(epsilon_weighted);
    k_NbN = n_NbN*k;
    k_NbO = n_NbO*k;
    
    for count_thickness=0:1:max_thickness
        d_NbN = count_thickness/5;
    
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
        for pass=1:1:100%number of internal reflections; generally three is adequate, but the calculation is fast
    
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
            omega = 2*pi*c/(lambda0*10^-9);
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
        T_plot(count_thickness+1, count_fill+1)=T_total;
        R_plot(count_thickness+1, count_fill+1) = R_total;
        A_plot(count_thickness+1, count_fill+1) = A_total;
        A_NbN_plot(count_thickness+1, count_fill+1) = A_NbN_total;
        thickness(count_thickness+1)=d_NbN;
        fill(count_fill+1)=fill_factor;
        R_plot+T_plot+A_plot; %check to make sure they add up to 100%
    end
end


%%
% surf(thickness, fill, A_NbN_plot*100, 'EdgeColor', 'None')
% view(2)
% colorbar
% caxis([0,100])
% xlabel('NbN thickness (nm)')
% ylabel('fill factor (%)')
% h=colorbar();
% ylabel(h, 'absorptance (%)')

surf(fill*100, thickness, A_NbN_plot*100, 'EdgeColor', 'None')
view(2)
colorbar
caxis([0,100])
xlabel('fill factor (%)')
ylabel('thickness (nm)')
h=colorbar();
ylabel(h, 'absorptance (%)')

thick2nm = A_NbN_plot(11,:);
thick4nm = A_NbN_plot(21,:);
thick6nm = A_NbN_plot(31,:);

figure
plot(fill*100, thick2nm*100, 'r-')
hold on
plot(fill*100, thick4nm*100, 'g.')
hold on
plot(fill*100, thick6nm*100, 'b--')
xlabel('fill factor (%)')
ylabel('absorptance (%)')