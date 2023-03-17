clear all
close all
clc

m_shortstack = readmatrix("02_12_2022_Run.csv");

test_time = m_shortstack(:,1);
start_time = test_time(1);
test_time = test_time - start_time;
voltage = m_shortstack(:,8);
xpos = [1:1:length(voltage)]';
subplot(3,1,1), plot(test_time, voltage, 'LineWidth',3.0,'Color',[0.4 0.2 0.6] )
hold on
grid on
xlabel('Time [s]');
ylabel('Voltage [V]');
grid minor
subplot(3,1,2), plot(xpos, voltage )
grid on
xlabel('Tabel idx [-]');
ylabel('Voltage [V]');
grid minor
hold on
switch_pos = [13 320 484 661 784 917 1013 1119 1195 1281 1337];
% xline(13,'LineStyle','--','LineWidth',2.0,'Alpha',0.5,'Color',[0.4 0.6 1.0])
% xline(320,'LineStyle','--','LineWidth',2.0,'Alpha',0.5,'Color',[0.4 0.6 1.0])
% xline(484,'LineStyle','--','LineWidth',2.0,'Alpha',0.5,'Color',[0.4 0.6 1.0])
% xline(661,'LineStyle','--','LineWidth',2.0,'Alpha',0.5,'Color',[0.4 0.6 1.0])
% xline(784,'LineStyle','--','LineWidth',2.0,'Alpha',0.5,'Color',[0.4 0.6 1.0])
% xline(917,'LineStyle','--','LineWidth',2.0,'Alpha',0.5,'Color',[0.4 0.6 1.0])
% xline(1013,'LineStyle','--','LineWidth',2.0,'Alpha',0.5,'Color',[0.4 0.6 1.0])
% xline(1119,'LineStyle','--','LineWidth',2.0,'Alpha',0.5,'Color',[0.4 0.6 1.0])
% xline(1195,'LineStyle','--','LineWidth',2.0,'Alpha',0.5,'Color',[0.4 0.6 1.0])
% xline(1281,'LineStyle','--','LineWidth',2.0,'Alpha',0.5,'Color',[0.4 0.6 1.0])
% xline(1337,'LineStyle','--','LineWidth',2.0,'Alpha',0.5,'Color',[0.4 0.6 1.0])

xline(switch_pos,'LineStyle','--','LineWidth',2.0,'Alpha',0.5,'Color',[0.4 0.6 1.0])

charge_current = 45.6;
discharge_current = -charge_current;
current_1 = ones(1,(switch_pos(1) ))*0.0;
current_2 = ones(1,-switch_pos(1)+ (switch_pos(2) ))*charge_current;
current_3 = ones(1,-switch_pos(2)+ (switch_pos(3) ))*discharge_current;
current_4 = ones(1,-switch_pos(3)+ (switch_pos(4) ))*charge_current;
current_5 = ones(1,-switch_pos(4)+ (switch_pos(5) ))*discharge_current;
current_6 = ones(1,-switch_pos(5)+ (switch_pos(6) ))*charge_current;
current_7 = ones(1,-switch_pos(6)+ (switch_pos(7) ))*discharge_current;
current_8 = ones(1,-switch_pos(7)+ (switch_pos(8) ))*charge_current;
current_9 = ones(1,-switch_pos(8)+ (switch_pos(9) ))*discharge_current;
current_10 = ones(1,-switch_pos(9)+ (switch_pos(10) ))*charge_current;
current_11 = ones(1,-switch_pos(10)+ (switch_pos(11) +1))*discharge_current;
combined_current = [current_1 current_2 current_3 current_4 current_5 current_6 current_7 current_8 current_9 current_10 current_11]';

subplot(3,1,3), plot(xpos, combined_current )
grid on
xlabel('Tabel idx [-]');
ylabel('Current [A]');
xline(switch_pos,'LineStyle','--','LineWidth',2.0,'Alpha',0.5,'Color',[0.4 0.6 1.0])

% time_step = 4;
% from_idx = floor(4 / time_step);
% to_idx = floor( 1936 / time_step);
% 
% one_cycle_time = test_time(from_idx:to_idx);
% one_cycle_voltage = voltage(from_idx:to_idx);
% 
% figure
% plot(one_cycle_time, one_cycle_voltage);
% grid on
% 

figure
time_step = 4;
max_time = test_time(length(test_time))
min_time = test_time(1)
step_time = (max_time - min_time) / length(test_time)
stop_idx = switch_pos(3); %one cycle

test_time_csv = test_time(1:stop_idx);
voltage_csv = voltage(1:stop_idx);
current_csv = combined_current(1:stop_idx);
subplot(2,1,1), plot(test_time_csv, voltage_csv, 'LineWidth',3.0,'Color',[0.4 0.2 0.6] )
hold on
grid on
xlabel('Time [s]');
ylabel('Voltage [V]');
grid minor

subplot(2,1,2), plot(test_time_csv, current_csv )
grid on
xlabel('Time [s]');
ylabel('Current [A]');
grid minor

system('echo time,voltage,current> batterydata.csv')
newcsv = [test_time_csv' ; voltage_csv'; current_csv']';
writematrix(newcsv,'batterydata.csv','WriteMode','append');


