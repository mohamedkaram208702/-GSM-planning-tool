function [Cluster_Size, Cells_Number, Cell_radius, T_intensity_cell, T_intensity_sector, Power_TransBS] = Network_Planning(GOS, City_Area, User_Density, SIR_min, sectorization_angle, flag)
%% Constants

n = 4;                                      % Loss exponant
c = 3*1e8;                                  % Speed of Light
f = 900e6;                                  % Frequency
hBS = 20;                                   % Base station height
hMS = 1.5;                                  % Mobile stattion height
PRX_min_dBm = -95;                          % Receiver sensitivity in dBm
S = 340;                                    % no. of channels
io = (sectorization_angle/360)*6;           % no. of interference cells
Sectors_number = 360/sectorization_angle;   % no. of sectors per cell
Au = 0.025;                                 % User traffic intensity

%% System Design

N = ((((10^(SIR_min/10))*io))^(2/n))/3;     % no. of cells per cluster

[i,k] = meshgrid(0:1:100,0:1:100);
N1 = i.^2 + i.*k + k.^2;
N1 = reshape(N1, 1, []);
N1 = unique(sort(N1));
N1 = N1(2:end);

for x = 1: length(N1)
    if N1(x)>= N
        N = N1(x);
        break;
    end
end
Cluster_Size = N;

if (N==1)&&(sectorization_angle == 360)
    Cells_Number=0;
    Cell_radius=0;
    T_intensity_cell=0;
    T_intensity_sector=0;
    Power_TransBS=0;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C_sector = floor(S/(N*Sectors_number));
k = Sectors_number*C_sector;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms x K
Summation = symsum((x.^K)/factorial(K),K,0,C_sector);
solx = vpasolve((GOS/100)==((x.^C_sector)/(Summation*factorial(C_sector))),x);
A_sector=solx(real(solx)>0&imag(solx)==0);
T_intensity_sector = A_sector;
T_intensity_cell = A_sector*Sectors_number;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_sector = floor(A_sector/Au);
area_sector = u_sector/User_Density;           % Area of sector, UD is the user density
Sectors_number = 360/sectorization_angle;       %n_sector = no. of sectors per cell, sectorization_angle is an input, 
area_cell = Sectors_number*area_sector;        % Area of one cell
R = sqrt((2/(3*sqrt(3))*area_cell));
Cell_radius = eval(R);
area_cell = Sectors_number*area_sector;        % Area of one cell
MN = ceil(City_Area/area_cell);                % MN is the total no. of cells, City_Area is an input
M = eval(MN/N);                                      % M is the number of clusters

Cells_Number = eval(MN);

CH = 0.8 + (1.1*log10(f*1e-6) - 0.7)*hMS - 1.56*log10(f*1e-6);                                      % Correction factor of height
Lv = 69.55 + 26.16*log10(f*1e-6) - 13.82*log10(hBS) - CH + (44.9 - 6.55*log10(hBS))*log10(R*1e-3);  % Losses

PTX = PRX_min_dBm + abs(Lv);
Power_TransBS = eval(PTX);


if flag == 1               % This flag is used to neglect power plots during verification when 
                           % "Network_Planning" function is called
%% PLOTS
d = 0:1:10e3;
PRX = PTX - (69.55 + 26.16*log10(f*1e-6) - 13.82*log10(hBS) - CH + (44.9 - 6.55*log10(hBS))*log10(d*1e-3));
plot (d, PRX ,'linewidth',1.5);
set(gca, 'XScale', 'log')
title('MS Received Power VS Distance from the BS using Hata Model','linewidth',1.5,'FontName','Times')
xlabel('Receiver distance from the BS (m)','FontName','Times')
ylabel('MS Received Power (dBm)','FontName','Times')
grid on

d_break = (4*hBS*hMS)*f/(c);
d = 0:1:d_break;

P_Friis = PTX + 20*log10(c./(4*pi*f*d));
figure;
plot(d,P_Friis,'black','linewidth',1.5);
set(gca, 'XScale', 'log')
grid on
hold on

d = d_break:1:10e3;
P_breakpoint = PTX + 20*log10(c./(4*pi*f*d_break)) - 10*n*log10(d/d_break);
plot(d,P_breakpoint,'blue','linewidth',1.5);
hold on
set(gca, 'XScale', 'log')
text(d_break+50,-110,'(d_{break})','FontName','Times')


title('MS Received Power VS Distance (Breakpoint Model and 2-Ray Model)','FontName','Times')
xlabel('Receiver distance from the BS (m)','FontName','Times')
ylabel('MS Received Power (dBm)','FontName','Times')


P_2_ray = 10*log10(PTX*((hBS*hMS./d.^2)).^2);
plot(d,P_2_ray, 'red','linewidth',1.5);
hold on
xline(d_break,'--');
legend('Friis'' Law', 'Breakpoint Model', '2-Ray Model')
end
end