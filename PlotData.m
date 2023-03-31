clear all
close all
load('WouterDataFile.mat')

%%

startval = 1;
% stopval = 48*3600;
stopval = length(foo.TestTime);

I = foo.Current_mA(startval:stopval);
I = rmmissing(I);
coulcount = cumsum(I)/3600;
% ti = hours(seconds(1:length(I)));
ti = 1:length(I);

V = foo.Voltage_V(startval:stopval);
V = rmmissing(V);
V = rmoutliers(V);
V(V <= 0.3) = [];
V(V >= 1) = [];
% tV = hours(seconds(1:length(V)));


tV = 1:length(V);


cap = foo.Capacity_mAh(startval:stopval);
cap = rmmissing(cap);
cap = rmoutliers(cap);
maxvals = findpeaks(cap);

t = foo.Record(1:stopval);
%%

close all
figure(1)
subplot(2,1,1)
plot(tV,V)
hold on
% yline(0.659+0.02*0.659,'--r')
% yline(0.659-0.02*0.659,'--r')
% xline(25962,'--k')
% xline(25969,'--k')
ylim([0.2,1])
% xlim([25955 26000])
title('Response time')
ylabel('Cell voltage [V]')
subplot(2,1,2)
plot(ti,I)
% ylim([-1.1*max(I)-5,1.1*max(I)])
% xlim([25955 26000])
xlabel('Time [s]')
ylabel('Current [mA]')

% figure(2)
% plot(1:length(cap),cap)
% xlim([0 length(cap)])
% figure(3)
% plot(1:length(coulcount),coulcount)
% xlim([0 length(coulcount)])
% 
% figure(5)
% subplot(2,1,1)
% stem(1:36,maxvals(1:2:end))
% xlim([1 36])
% title('Peak cycle capacities')
% ylabel('Capacity [mAh]')
% subplot(2,1,2)
% stem(1:36,maxvals(1:2:end)/max(maxvals))
% ylabel('Normalized capacity')
% xlabel('Cycle count')
% xlim([1 36])
