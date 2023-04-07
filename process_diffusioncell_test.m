clear all
close all
clc

% If loading old data
load('WouterDataFile.mat')

% If loading new data
% foo = readtable('newOGmem.csv');

startval = 1;
stopval = length(foo.TestTime);

Ts = 1;

current = foo.Current_mA(startval:stopval)/1000;
current = rmmissing(current);
test_time = 1:length(current);
test_time = test_time';

voltage = foo.Voltage_V(startval:stopval);
voltage = rmmissing(voltage);

for ii = 1:length(voltage)
    if voltage(ii) <= 0.1
        voltage(ii) = voltage(ii-1);
    elseif voltage(ii) >= 1
        voltage(ii) = voltage(ii-1);
    else
    end
end
% voltage = rmoutliers(voltage);
% voltage(voltage <= 0.3) = [];
% voltage(voltage >= 1) = [];


% [current, TimeRegular] = resample(current,test_time,1/Ts,'linear');
% [voltage, test_time] = resample(voltage,test_time,1/Ts,'linear');


dsfac = 120;

test_time = downsample(test_time,dsfac);
voltage = downsample(voltage,dsfac);
current = downsample(current,dsfac);

figure()
nexttile
plot(test_time,voltage)
nexttile
plot(test_time,current)
%%

start_idx = 1;
% start_idx = 2480;
% start_idx = 2330;
% stop_idx = length(test_time);
% stop_idx = 35000/Ts;
stop_idx = 4*24*3600/dsfac;
% stop_idx = 510;
if stop_idx > length(voltage)
    stop_idx = length(voltage);
end
% stop_idx = 7100;
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

subplot(2,1,2), plot(test_time_csv, current_csv)
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
