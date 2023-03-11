close all
clear all
clc

%read simulationdata from LTSpice and make regular sampled interval. Add
%noise to "simulate" measurement noise
%m=csvread('RC_load_first_order.txt', 2);
m2 = readmatrix("RC_load_first_order_1_F_10_Ohm.txt");
m = rmmissing(m2);
%m = m{:,:}; %convert to matrix
voltage = m(:,2);
current = m(:,3);
irregTx = m(:,1);
desiredTs = 1;
desiredFs = 1/desiredTs ;
[voltage_regular, TimeRegular] = resample(voltage,irregTx,desiredFs,'linear');
[current_regular, TimeRegular] = resample(current,irregTx,desiredFs,'linear');

r1 = rand(1000,1);

voltage_noise = 0.2; %amplitude relative to mean
current_noise = 0.002; %amplitude relative to mean
v_n = rand(length(voltage_regular),1)*(2*voltage_noise) - voltage_noise;
i_n = rand(length(current_regular),1)*(2*current_noise) - current_noise;
voltage_noised =voltage_regular + v_n;
current_noised =current_regular + i_n;

subplot(2,1,1), plot(irregTx,voltage,'.-',TimeRegular,voltage_noised,'-')
legend('Original','Resampled')
grid on
subplot(2,1,2), plot(irregTx,current,'.-',TimeRegular,current_noised,'-')
legend('Original','Resampled')
grid on
newcsv = [TimeRegular' ; voltage_noised'; current_noised']';
system('echo time,voltage,current> batterydata.csv')
writematrix(newcsv,'batterydata.csv','WriteMode','append');
%csvwrite('batterydata.csv',newcsv);


%runs Rust program with 'batterydata.csv' -> takes time!!!
[status, result]= system('cargo run')

%post Genetic algorithm plots
m3 = readmatrix("fitness.csv");
figure
plot(m3(:,1), m3(:,2),'LineWidth',3,'Color',[0.6 0.2 0.3])
grid on
xlabel('Generations')
ylabel('Fitness')

simdata = readmatrix('batterydata.csv');
ga_data = readmatrix('best_result.csv');
simtime = simdata(:,1);
simvoltage = simdata(:,2);

gatime = ga_data(:,1);
gavoltage = ga_data(:,3);
figure
plot(simtime,simvoltage,'.-',gatime,gavoltage,'-')
legend('Input','GA estimated result')
grid on
xlabel('Time [s]')
ylabel('Voltage [V]')
