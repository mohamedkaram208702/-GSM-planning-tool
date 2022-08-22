flag = 0;               % This flag is used to neglect power plots during verification when 
                        % "Network_Planning" function is called
                        
Network_Planning(2, 100e6, 1400e-6, 19, 120,1);

%% part 1

% default numbers changed by calling the function if needed

City_Area=100*10^6;
User_Density=1400*10^-6;
SIR_min= 1:1:30;
sectorization_angle = [360 120 60];
N=zeros(size(SIR_min));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
for i = 1:3
for x = 1:length(SIR_min)
[Cluster_Size, Cells_Number, Cell_radius, T_intensity_cell, T_intensity_sector, Power_TransBS] = Network_Planning(2, 100e6, 1400e-6, SIR_min(x), sectorization_angle(i),flag);
N(x)=Cluster_Size;
end
hold on
grid on
plot(SIR_min,N,"linewidth",1.5);
title("Cluster size vs SIR_m_i_n @ Different Sectorization angles");
xlabel("SIR_m_i_n (dB)");
ylabel("Cluster size (N)");
end
legend('360^o sectorization','120^o sectorization','60^o sectorization');
hold off

%% parts 2 & 3 

GOS_2and3=1:1:30;
Cells_number=zeros(size(GOS_2and3));
A_cell=zeros(size(GOS_2and3));
for SIR_2and3 = [19 14]
figure;
for i = 1:3
for x = 1:length(GOS_2and3)
[Cluster_Size, Cells_Number, Cell_radius, T_intensity_cell, T_intensity_sector, Power_TransBS] = Network_Planning(GOS_2and3(x), 100*10^6, 1400e-6, SIR_2and3, sectorization_angle(i),flag);
Cells_number(x) = Cells_Number;
A_cell(x) = T_intensity_cell;
end

subplot(2,1,1);
hold on
grid on
plot(GOS_2and3,Cells_number,"linewidth",1.5);
title("Number of cells vs GOS @ SIR_m_i_n = "+ SIR_2and3);
xlabel("GOS (%)");
ylabel("Number of cells");
if i == 3
    legend('360^o sectorization','120^o sectorization','60^o sectorization');
end

subplot(2,1,2);
hold on
grid on
plot(GOS_2and3,A_cell,"linewidth",1.5);
title("Traffic Intensity per cell vs GOS @ SIR_m_i_n = "+ SIR_2and3);
xlabel("GOS (%)");
ylabel("Traffic Intensity per cell");
if i == 3
    legend('360^o sectorization','120^o sectorization','60^o sectorization');
end
end
hold off
end

%% 4 & 5

angle = [360 120 60];
SIR = [14 19];
Density = 100:100:2000;
N = zeros(size(Density));
R = zeros(size(Density));

for i = 1:1:2
    figure;
    grid on
    SIR_i = SIR(i);
        for j = 1:1:3
            sectorization_angle = angle(j);
            for k = 1:length(Density)
                [Cluster_Size, Cells_Number, Cell_radius, T_intensity_cell, T_intensity_sector, Power_TransBS] = Network_Planning(2, 100e6, Density(k)*1e-6, SIR_i, sectorization_angle, flag);
                N(k) = Cells_Number;
                R(k) = Cell_radius;
            end
            hold on
            plot(Density, N, "linewidth",1.5)
            legend("Sectorization at " + angle(1) + " degrees", "Sectorization at " + angle(2) + " degrees", "Sectorization at " + angle(3) + " degrees")
        end
            title("Total Number of Cells in System VS User Density at SIR = " + SIR(i) + " dB",'FontName','Times')
            xlabel('User Density (users/km^{2})','FontName','Times')
            ylabel('Total Number of Cells in System','FontName','Times')
end


for i = 1:1:2
    SIR_i = SIR(i);
    figure;
    grid on

        for j = 1:1:3
            sectorization_angle = angle(j);
            for k = 1:length(Density)
                [Cluster_Size, Cells_Number, Cell_radius, T_intensity_cell, T_intensity_sector, Power_TransBS] = Network_Planning(2, 100e6, Density(k)*1e-6, SIR_i, sectorization_angle,flag);
                N(k) = Cells_Number;
                R(k) = Cell_radius;
            end
            hold on
            plot(Density, R, "linewidth",1.5)
            legend("Sectorization at " + angle(1) + " degrees", "Sectorization at " + angle(2) + " degrees", "Sectorization at " + angle(3) + " degrees")
        end
            title("Cell Radius VS User Density at SIR = " + SIR(i) + " dB",'FontName','Times')
            xlabel('User Density (users/km^{2})','FontName','Times')
            ylabel('Cell Radius','FontName','Times')
end