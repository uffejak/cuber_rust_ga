
%runs Rust program with 'batterydata.csv' -> takes time!!!
%[status, result]= system('cargo run')
close all
clear all
clc


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
gacharge = ga_data(:,2);
gavoltage = ga_data(:,3);
figure
subplot(2,1,1), plot(simtime,simvoltage,'.-',gatime,gavoltage,'-')
legend('Input','GA estimated result')
grid on
xlabel('Time [s]')
ylabel('Voltage [V]')
hold on
subplot(2,1,2), plot(gatime,gacharge,'-')
legend('GA estimated result')
grid on
xlabel('Time [s]')
ylabel('Charge [C]')
hold on
