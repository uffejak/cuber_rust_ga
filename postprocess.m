
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
gac1a = ga_data(:,4);
gac1c = ga_data(:,5);
gac2a = ga_data(:,6);

gac0c = 2000-gac1a-gac1c-gac2a;
gatime = ga_data(:,1);
gacharge = ga_data(:,2);
gavoltage = ga_data(:,3);
gacurrent = ga_data(:,end);

figure(1)
subplot(2,1,1)
plot(simtime,simvoltage,'.-',gatime,gavoltage,'-')
legend('Measured Voltage','Estimated Voltage')
grid on
xlabel('Time [s]')
ylabel('Voltage [V]')
xlim([gatime(1) gatime(end)])
subplot(2,1,2)
plot(gatime,gacurrent)
grid on
xlabel('Time [s]')
ylabel('Current [A]')
xlim([gatime(1) gatime(end)])

figure(2)
plot(gatime,gac1a,'--r')
hold on
plot(gatime,gac2a,'--k')
plot(gatime,gac1c,'--b')
hold off
xlim([gatime(1) gatime(end)])
title('Tank Concentrations')
xlabel('Time [s]')
ylabel('Concentration [$\frac{mol}{m^3}$]')
legend('$Cu^{1+}$ (anolyte)','$Cu^{2+}$ (anolyte)','$Cu^{1+}$ (catholyte)') 

figure(3)
plot(gatime,gac1a+gac2a+gac1c+gac0c)