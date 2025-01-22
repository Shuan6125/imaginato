clear all
clc

%% Define the dimensions of the timber
depth = 0.025; % Depth of the timber in meters

% Define the properties of the timber
E12 = 10.3e9*0.154; % radius MOE of oak at MC=12 (Pa)
Eg = 7.9e9*0.154; %radius MOE of oak at green status (Pa)
MC_0 = 35; % Initial moisture content(%) as a fraction of the oven dry weight

% Define the boundary conditions
Temp=[37.78 43.3 46.1 48.89 54.4 60 65.56];
RH=[80 73 66 58 42 33 30];
EMC=FindEMC(Temp,RH); %EMC for each stage
Time=[2 4 2 2 2 2 2]; %time period in days
    if length(EMC)~=length(Time)
        error("schedule length not equal")
    end

% Define the grid for the moisture content
NGrid = 251; %Number of spatial points, = NumEle+1
z = linspace(0, depth, NGrid); % Spatial grid
dz= abs(z(2)-z(1));
MC = zeros(1, length(z)); % Moisture content
MC(:) = MC_0; % Initialize the moisture content

% Define the time step and number of time steps
dt = 1; % Time step in seconds
tDay = 86400/dt; % Num of time step per day, 1 Day = 86400 s
NRec = 240; %Number of data recording per day, time step number for stress calculation
NUp = 24; %Number of Du updating per day
tRec = round(tDay/NRec); % time step interval for data recording
tUp=round(tDay/NUp);% time step interval for Du updating.

TotTime=sum(Time);
nt=TotTime*tDay;

%Initializing Record
RecordMC=zeros(floor(nt/tRec),length(z));
tRecord=0;
nRecord=1;
RecordDu=zeros(floor(nt/tUp),length(z));
nDu=1;

LSch=length(EMC);
for i=1:LSch

BC_EMC = EMC(i); % Ambient EMC in degrees Celsius
N = Time(i)*tDay; % Number of time steps for this stage
MC(end) = BC_EMC;
MC(1) = BC_EMC ;
D_u=19.2*exp(-5280/(273+Temp(i)))*exp(1.15*MC./100)*0.0001;

for n = 1:N

dU = D_u(2:end-1) .* dt .* (MC(3:end) - 2 * MC(2:end-1) + MC(1:end-2)) ./ (dz^2);
MC(2:end-1) = MC(2:end-1) + dU;
MC(end) = BC_EMC;
MC(1) = BC_EMC ;
    if mod(n,tRec)==0
        RecordMC(nRecord,:)=MC;
        tRecord=tRecord+tRec;
        nRecord=nRecord+1;
    end
    if mod(n,tUp)==0
    RecordDu(nDu,:)=D_u;
    D_u=19.2*exp(-5280/(273+Temp(i)))*exp(1.15*MC./100)*0.0001;% Oak desorption coefficient 
    nDu=nDu+1;
    end
end

end

%% Plot the moisture content along the thickness of the timber
figure;
plot(z, MC);
ylabel('Moisture content (%)');
xlabel('Position along the thickness of the lumber (m)');
title('Final Moisture Content Distribution');

[m,n]=size(RecordMC);
figure;
hold on;
for e=1:60:m
    plot(z,RecordMC(e,:));
end
ylabel('Moisture content (%)');
xlabel('Position along the thickness of the lumber (m)');
title('Moisture Content Profiles Over Time');

% Moisture content surface plot
figure;
TimeVector = (0:tRec:(size(RecordMC, 1)-1)*tRec) / 86400; % Time in days
[X, Y] = meshgrid(z, TimeVector);
surf(X, Y, RecordMC, 'EdgeColor', 'none');
shading interp;

% Add custom grid lines for the first surf plot
hold on;
position_step = round(length(z)/ 10); % Position grid step (10 divisions)
time_step = round((24 * 3600) / tRec); % Time grid step (24 hours)

% Add position grid lines
for i = 1:position_step:length(z)
    plot3(z(i)*ones(size(TimeVector)), TimeVector, RecordMC(:, i), 'k');
end

% Add time grid lines
for j = 1:time_step:size(TimeVector, 2)
    plot3(z, TimeVector(j)*ones(size(z)), RecordMC(j, :), 'k');
end

xlabel('Position along thickness (m)');
ylabel('Time (days)');
zlabel('Moisture content (%)');
title('Moisture Content Over Time (Surface Plot with Custom Grids)');
colorbar;

%% Calculate drying stress
MCTimeSeries=transpose(RecordMC);
alpha = 0.0011; % Coefficient of hygroscopic expansion (strain per unit moisture content, unit: %^-1)
m= 0.085;
b= -0.045; %Mpa^-1

[N, T] = size(MCTimeSeries);

% Initialize the drying stress matrix
drying_strain_HE = zeros (N,T+1);
Ele_drying_strain_HE = zeros (N-1,T+1);
Factor_MS = zeros (N,T+1);
Ele_Factor_MS = zeros (N-1,T+1);
drying_stress = zeros(N-1, T+1);
drying_stress_inc = zeros(N-1,T+1);
Em=zeros(T+1,N);
Ele_Em=ones(T+1,N-1)*350e6;
moisture_gradient = diff(MCTimeSeries, 1, 1) / dz;
moisture_change=transpose(diff(transpose(MCTimeSeries)));
moisture_change=[zeros(N,1) moisture_change];
x=(z(2:end)+z(1:end-1))/2;
RecSol=zeros(2,T+1);
RecEq11=zeros(1,T+1);
RecEq12=zeros(1,T+1);
RecEq13=zeros(1,T+1);
RecEq21=zeros(1,T+1);
RecEq22=zeros(1,T+1);
RecEq23=zeros(1,T+1);

for t = 2:T+1
    drying_strain_HE(:, t) =  alpha * moisture_change(:, t-1);
    Ele_drying_strain_HE = (drying_strain_HE(1:end-1, t)+drying_strain_HE(2:end, t))/2;
    Factor_MS(:, t) = m*abs(moisture_change(:, t-1)/100)+b*moisture_change(:, t-1)/100;
    Ele_Factor_MS(:, t) = (Factor_MS(1:end-1,t)+Factor_MS(2:end,t))/2;

    eq11=sum(dz.*Ele_Em(t,:));
    eq12=sum((x-0.025/2).*dz.*Ele_Em(t,:));
    eq13=sum(dz.*Ele_Em(t,:).*(Ele_drying_strain_HE+drying_stress(:, t-1)*(10^(-6)).*Ele_Factor_MS(:, t))');

    eq21=sum(x.*dz.*Ele_Em(t,:));
    eq22=sum(x.*(x-0.025/2).*dz.*Ele_Em(t,:));
    eq23=sum(x.*dz.*Ele_Em(t,:).*(Ele_drying_strain_HE+drying_stress(:, t-1)*(10^(-6)).*Ele_Factor_MS(:, t))');

    A=[eq11 eq12; eq21 eq22];
    B=[eq13; eq23];
    Sol=linsolve(A,B);
    RecSol(:,t)=Sol;

    drying_stress_inc(:, t) = Ele_Em(t,:)'.*(Sol(1,1)+Sol(2,1)*(transpose(x-0.025/2))-Ele_drying_strain_HE-drying_stress(:, t-1)*(10^(-6)).*Ele_Factor_MS(:, t));
    drying_stress(:, t)=drying_stress(:, t-1)+drying_stress_inc(:, t);
end

DryingStress=transpose(drying_stress);
DryingStressR=DryingStress(2:end,:);
figure;
hold on;
imagesc(DryingStress);
colorbar;
ylabel('Time step');
xlabel('Spatial position (element)');
title('Drying Stress Distribution (Pa)');

DryingStressGPa=DryingStress/1000000;
figure;
surfc(DryingStressGPa, 'EdgeColor', 'none');
shading interp;

% Add custom grid lines for the second surf plot
hold on;
position_step = 25; % Position grid step (10 divisions)
time_step = round((24 * 3600) / tRec); % Time grid step (24 hours)

% Add position grid lines
for i = 1:position_step:250
    plot3(x(i)*ones(1, size(DryingStressGPa, 1))*10000, 1:size(DryingStressGPa, 1), DryingStressGPa(:, i), 'k');
end

% Add time grid lines
for j = 1:time_step:size(DryingStressGPa, 1)
    plot3(x*10000, j*ones(1, length(x)), DryingStressGPa(j, :), 'k');
end

colorbar;
ylabel('Time step');
xlabel('Thickness position (0.1 mm)');
title('Drying Stress (MPa) with Custom Grids');

figure;
plot(transpose(x),DryingStress(end,:));
xlabel('Position along thickness (m)');
ylabel('Stress (Pa)');
title('Final Drying Stress Distribution');

MaxDryingStress=zeros(1,NGrid-1);
MinDryingStress=zeros(1,NGrid-1);
for i=1:NGrid-1
     MaxDryingStress(1,i)=max(drying_stress(i,:));
     MinDryingStress(1,i)=min(drying_stress(i,:));   
end
figure;
plot(x,MaxDryingStress, '-r', 'DisplayName', 'Max Stress');
hold on;
plot(x,MinDryingStress, '-b', 'DisplayName', 'Min Stress');
xlabel('Position along thickness (m)');
ylabel('Stress (Pa)');
legend;
title('Max and Min Stress Distribution');

tRecord = (0:tRec:nt)/86400;
TimeMaxDryingStress=zeros(1,length(tRecord));
TimeMinDryingStress=zeros(1,length(tRecord));
for t=1:length(tRecord)
    TimeMaxDryingStress(t)=max(drying_stress(:,t));
    TimeMinDryingStress(t)=min(drying_stress(:,t));
end
figure;
plot(tRecord,TimeMaxDryingStress, '-r', 'DisplayName', 'Max Stress');
hold on;
plot(tRecord,TimeMinDryingStress, '-b', 'DisplayName', 'Min Stress');
xlabel('Time (days)');
ylabel('Stress (Pa)');
legend;
title('Max and Min Stress Over Time');
