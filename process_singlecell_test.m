clear all
close all
clc

m_shortstack = readmatrix("WouterCellData.csv");
% m_shortstack = readmatrix("WouterCellData.xlsx");

test_time = m_shortstack(:,3);
start_time = test_time(1);
test_time = test_time - start_time;
voltage = m_shortstack(:,1);
current = m_shortstack(:,2);


dsfac = 10/5;

current = downsample(current,dsfac);
voltage = downsample(voltage,dsfac);
test_time = 0:5*dsfac:length(voltage)*5*dsfac-5*dsfac;
test_time = test_time';

% [voltage, TimeRegular] = resample(voltage,test_time,1/Ts,'linear');
% [current, test_time] = resample(current,test_time,1/Ts,'linear');

%%
%voltage_offset = voltage(1);
%voltage = voltage - voltage_offset;
%flow = m_shortstack();
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
xline(switch_pos,'LineStyle','-.' ,'LineWidth',2.0,'Alpha',0.5,'Color',[0.4 0.6 1.0])

chargestart_pos   = [023 305 490 650 790 910 1020 1110 1205 1270];
dischargestart_pos= [328 475 667 770 921 995 1125 1180 1285 1325];
xline(chargestart_pos,'LineStyle',':','LineWidth',2.0,'Alpha',0.5,'Color',[1.0 0.4 0.6])
xline(dischargestart_pos,'LineStyle','--','LineWidth',2.0,'Alpha',0.5,'Color', [0.2 0.7 0.3])

% charge_current = 45.6;
% discharge_current = -charge_current*1.1;
% current_1 = ones(1,(switch_pos(1) ))*0.0;
% current_2 = ones(1,-switch_pos(1)+ (switch_pos(2) ))*charge_current;
% current_3 = ones(1,-switch_pos(2)+ (switch_pos(3) ))*discharge_current;
% current_4 = ones(1,-switch_pos(3)+ (switch_pos(4) ))*charge_current;
% current_5 = ones(1,-switch_pos(4)+ (switch_pos(5) ))*discharge_current;
% current_6 = ones(1,-switch_pos(5)+ (switch_pos(6) ))*charge_current;
% current_7 = ones(1,-switch_pos(6)+ (switch_pos(7) ))*discharge_current;
% current_8 = ones(1,-switch_pos(7)+ (switch_pos(8) ))*charge_current;
% current_9 = ones(1,-switch_pos(8)+ (switch_pos(9) ))*discharge_current;
% current_10 = ones(1,-switch_pos(9)+ (switch_pos(10) ))*charge_current;
% current_11 = ones(1,-switch_pos(10)+ (switch_pos(11) +1))*discharge_current;
% combined_current = [current_1 current_2 current_3 current_4 current_5 current_6 current_7 current_8 current_9 current_10 current_11]';

subplot(3,1,3), plot(xpos, current )
grid on
xlabel('Tabel idx [-]');
ylabel('Current [A]');
xline(switch_pos,'LineStyle','--','LineWidth',2.0,'Alpha',0.5,'Color',[0.4 0.6 1.0])

generate_sections = false;
figure
max_time = test_time(length(test_time));
min_time = test_time(1);
step_time = (max_time - min_time) / length(test_time);

% if generate_sections
%     start_idx = chargestart_pos(1);
%     stop_idx = chargestart_pos(2); %one charge cycle
% else
% %     start_idx = chargestart_pos(3);
% %     stop_idx = switch_pos(11); %one cycle
%     start_idx = chargestart_pos(1);
%     stop_idx = switch_pos(3); %one cycle
% end

% desiredTs = 4;
% desiredFs = 1/desiredTs ;
% [voltage_regular, TimeRegular] = resample(voltage,test_time,desiredFs,'linear');
% [current_regular, TimeRegular] = resample(combined_current,test_time,desiredFs,'linear');

start_idx = 1;
% start_idx = 9600/10;
% start_idx = 12680/10;
% start_idx = 2480;
% stop_idx = length(test_time);
% stop_idx = 35000/Ts;
% stop_idx = 4830;
stop_idx = length(voltage);
% stop_idx = 12660/10;
if stop_idx > length(voltage)
    stop_idx = length(voltage);
end
test_time_csv = test_time(start_idx:stop_idx);
voltage_csv = voltage(start_idx:stop_idx);
current_csv = current(start_idx:stop_idx);

subplot(2,1,1), plot(test_time_csv, voltage_csv)
hold on
xlim("tight")
xrange_used = xlim;
xlim(xrange_used);
xlim('manual'); %locks limit to current values
grid on
xlabel('Time [s]');
ylabel('Voltage [V]');
grid minor
%subplot(2,1,1), plot(test_time, voltage, 'LineWidth',1.0,'Color',[1.0 0.4 0.2],'LineStyle','--' )

subplot(2,1,2) 
plot(test_time_csv, current_csv)
hold on
grid on
xlabel('Time [s]');
ylabel('Current [A]');
grid minor
xlim("tight")
xlim(xrange_used);
xlim('manual'); %locks limit to current values
%subplot(2,1,2), plot(test_time, combined_current , 'LineWidth',1.0,'Color',[1.0 0.4 0.2],'LineStyle','--')

system('echo time,voltage,current> batterydata.csv')
newcsv = [test_time_csv' ; voltage_csv'; current_csv']';
    writematrix(newcsv,'batterydata.csv','WriteMode','append');

disp('Done');
