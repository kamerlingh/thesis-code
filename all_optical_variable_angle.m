%Multiple angle model based on Knittel (Optics of Thin Films)

%TM, parallel, p-polarization: electric field in plane of incidence
%TE, perpendicular, s-polarization: electric field perpendicular to plane
%of incidence
clear all
lambda0 = 1550; %all distances are in nanometers

n1 = 1.5;

n_NbN = 5.23-5.82*1i;
n_air =1;
fill_factor = 0.4;
epsilon_weighted_NbN = fill_factor*(n_NbN^2)+(1-fill_factor)*(n_air^2);
n2 = sqrt(epsilon_weighted_NbN);
d2 = 4;

imeps = (8.85418782*10^-12)*imag(n2^2);
c = 3*10^8; %c in m/s
omega = 2*pi*c/(lambda0*10^-9);

n3 = 2.217;
d3 = 20;
figure
for d3=0:10:160
    

n4 = 0.55-11.5*1i;
d4 = 120;
n5 = 1;

max = 90;

angle = zeros(max,1);
Rplot = zeros(max,1);
Tplot = zeros(max,1);
A_NbN = zeros(max,1);
A_gold = zeros(max,1);

%P-polarization math
for count = 1:1:max
    
   theta1 = (count-1)*pi/180;
   theta2 = asin((n1/n2)*sin(theta1));
   theta3 = asin((n1/n3)*sin(theta1));
   theta4 = asin((n1/n4)*sin(theta1));
   theta5 = asin((n1/n5)*sin(theta1));
   
   Y1 = n1/cos(theta1);
   Y2 = n2/cos(theta2);
   Y3 = n3/cos(theta3);
   Y4 = n4/cos(theta4);
   Y5 = n5/cos(theta5);
   
   phi2 = 2*pi*n2*d2*cos(theta2)/lambda0;
   phi3 = 2*pi*n3*d3*cos(theta3)/lambda0;
   phi4 = 2*pi*n4*d4*cos(theta4)/lambda0;
   
   W12 = 0.5*[1+(Y2/Y1) 1-(Y2/Y1); 1-(Y2/Y1) 1+(Y2/Y1)];
   U2 = [exp(1i*phi2) 0;0 exp(-1i*phi2)];
   W23 = 0.5*[1+(Y3/Y2) 1-(Y3/Y2); 1-(Y3/Y2) 1+(Y3/Y2)];
   U3 = [exp(1i*phi3) 0;0 exp(-1i*phi3)];
   W34 = 0.5*[1+(Y4/Y3) 1-(Y4/Y3); 1-(Y4/Y3) 1+(Y4/Y3)];
   U4 = [exp(1i*phi4) 0;0 exp(-1i*phi4)];
   W45 = 0.5*[1+(Y5/Y4) 1-(Y5/Y4); 1-(Y5/Y4) 1+(Y5/Y4)];
   
   M = W12*U2*W23*U3*W34*U4*W45;
   
   a=cos(theta1)/cos(theta5);
   
   R = abs(M(2,1)/M(1,1))^2;
   T = (n5*cos(theta5)/(n1*cos(theta1)))*abs(a/M(1,1))^2;
   
   %absorptance in NbN
   matrix = W23*U3*W34*U4*W45;
   AB = matrix*[abs(1/M(1,1)) 0]';
   d = d2*10^-9; %d is in meters
   X = 0:d/1000:d;
   k_NbNm = (n2*2*pi*cos(theta2))/(lambda0*10^9); %k_NbN in meters, not nm
   Y = abs(AB(1)*exp(-1i*k_NbNm*(X-d))+AB(2)*exp(1i*k_NbNm*(X-d))).^2;
   Z = trapz(X,Y);
   irradiance = n1*c*(8.85418782*10^-12)*0.5*cos(theta2)/cos(theta1);%need to include n_substrate for agreement with previous simulations
   Q = 0.5*omega*imeps*Z/irradiance;
   A_NbN(count) = -100*Q;
   
   angle(count) = theta1*180/pi;
   Rplot(count) = R*100;
   Tplot(count) = T*100;
   
end

%S-polarization math
% for count = 1:1:max
%     
%    theta1 = (count-1)*pi/180;
%    theta2 = asin((n1/n2)*sin(theta1));
%    theta3 = asin((n1/n3)*sin(theta1));
%    theta4 = asin((n1/n4)*sin(theta1));
%    theta5 = asin((n1/n5)*sin(theta1));
%    
%    Y1 = -n1*cos(theta1);
%    Y2 = -n2*cos(theta2);
%    Y3 = -n3*cos(theta3);
%    Y4 = -n4*cos(theta4);
%    Y5 = -n5*cos(theta5);
%    
%    phi2 = 2*pi*n2*d2*cos(theta2)/lambda0;
%    phi3 = 2*pi*n3*d3*cos(theta3)/lambda0;
%    phi4 = 2*pi*n4*d4*cos(theta4)/lambda0;
%    
%    W12 = 0.5*[1+(Y2/Y1) 1-(Y2/Y1); 1-(Y2/Y1) 1+(Y2/Y1)];
%    U2 = [exp(1i*phi2) 0;0 exp(-1i*phi2)];
%    W23 = 0.5*[1+(Y3/Y2) 1-(Y3/Y2); 1-(Y3/Y2) 1+(Y3/Y2)];
%    U3 = [exp(1i*phi3) 0;0 exp(-1i*phi3)];
%    W34 = 0.5*[1+(Y4/Y3) 1-(Y4/Y3); 1-(Y4/Y3) 1+(Y4/Y3)];
%    U4 = [exp(1i*phi4) 0;0 exp(-1i*phi4)];
%    W45 = 0.5*[1+(Y5/Y4) 1-(Y5/Y4); 1-(Y5/Y4) 1+(Y5/Y4)];
%    
%    M = W12*U2*W23*U3*W34*U4*W45;
%    
%    R = abs(M(2,1)/M(1,1))^2;
%    T = (n3*cos(theta3)/(n1*cos(theta1)))*abs(1/M(1,1))^2;
%    
%    %absorptance in NbN
%    matrix = W23*U3*W34*U4*W45;
%    AB = matrix*[abs(1/M(1,1)) 0]';
%    d = d2*10^-9; %d is in meters
%    X = 0:d/1000:d;
%    k_NbNm = (n2*2*pi*cos(theta2))/(lambda0*10^9); %k_NbN in meters, not nm
%    Y = abs(AB(1)*exp(-1i*k_NbNm*(X-d))+AB(2)*exp(1i*k_NbNm*(X-d))).^2;
%    Z = trapz(X,Y);
%    irradiance = n1*c*(8.85418782*10^-12)*0.5*cos(theta2)/cos(theta1);%need to include n_substrate for agreement with previous simulations
%    Q = 0.5*omega*imeps*Z/irradiance;
%    A_NbN(count) = -100*Q;
%    
%    angle(count) = theta1*180/pi;
%    Rplot(count) = R*100;
%    Tplot(count) = T*100;
   
% end

hold on
plot(angle, A_NbN, 'k-')
xlabel('angle (degrees)')
ylabel('absorptance in NbN (%)')
end