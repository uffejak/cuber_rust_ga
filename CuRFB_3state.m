clear all
close all
load('WouterCellData.csv')

%%

% Semi-believable values that make pretty graphs
k = 10^-12;
E0 = 0.65;
Ri = 0.2;
a = 0.5;
kp = 7.853*10^-7;
kn = 3.456*10^-6;

% Species ordering (soluble) is c1 c1 c2
D = [0 0 0; 0 0 2*k; 0 0 -k];

currentsigns = [-1;-1;1];

% Model parameters
V_tank = 50E-6; 
V_cell = 0.0025*0.002; 
N = 1; 
F = 96845; 
R = 8.314;  
z = 1;
S = 0.0025; 
d = 5*10^-5;
dt = WouterCellData(2,3)-WouterCellData(1,3);
molarity = 700;

c_t = [molarity;molarity;0];
c_c = [molarity;molarity;0];
c0 = 0;

Ival = 0.5;

NumSteps = 1*3600;
Q = 5e-7;
cntr = 0;
holup = 0;
waittime = 0;
I = WouterCellData(:,2);

% State dynamics of cells and tanks
fdyn = @(x,u) [x(1:3)+dt/V_cell*(u(1)*[1 0 0; 0 1 0; 0 0 1]*(x(4:6)-x(1:3))+1/(z*F)*currentsigns*u(2)+S/d*D*x(1:3));...
    x(4:6)+dt/V_tank*(N*u(1)*[1 0 0; 0 1 0; 0 0 1]*(x(1:3)-x(4:6)))];

% Exchange current functions

jpfun = @(x,u) 1/S*(F*kp*x(3)^(1-a)*x(1)^a);
jnfun = @(x,u) 1/S*(F*kn*x(2)^(a)*1000^(1-a));

% Voltages

Vpfun = @(x,u) 2*R*323.15/F*log(1/(2*jpfun(x,u)*S)*u(2)+sqrt((1/(2*jpfun(x,u)*S)*u(2))^2+1));
Vnfun = @(x,u) 2*R*323.15/F*log(1/(2*jnfun(x,u)*S)*u(2)+sqrt((1/(2*jnfun(x,u)*S)*u(2))^2+1));
Vnerfun = @(x,u) R*293.15/F*log((x(3)/x(1)+10e-12)*(1000/((x(2)+10e-12))));

% Full voltage dynamics
Vfun = @(x,u) N*(E0+Vnerfun(x,u)+Vpfun(x,u)-Vnfun(x,u)+Ri*u(2));


M = 100;
n = 6;
p = 1;
% myPF = SimplePF(M,n,p,fdyn,Vfun,[c_c;c_t]);
% xp(:,1) = [c_c;c_t];

for ii = 1:NumSteps
    SOC(ii) = c_t(3,ii)/(c_t(1,ii)+c_t(3,ii)); % SOC
    realSOC(ii) = c_t(3,ii)/(c_t(1,1)+c_t(3,1)); % SOC relative to initial
    
    % State dynamics
      states = [c_c(:,ii);c_t(:,ii)]; % Combine the states so they fit the function interface
      res = fdyn(states,[Q,I(ii)]); % Run function
      c_c(:,ii+1) = max(res(1:3),0);   % Split again for plotting
      c_t(:,ii+1) = max(res(4:6),0); 
    
    % Cu(0) from mass balance (in weird unit)
    c0(:,ii+1) = 2*molarity-sum(c_c(:,ii+1));

        if holup ~= 0
        I(ii+1) = 0;
        cntr = cntr + 1;
        if cntr > waittime
            I(ii+1) = sign(holup)*Ival;
        end
        if cntr > waittime+120
            holup = 0;
            cntr = 0;
        end
    else
        if (SOC(ii) > 0.99)
            holup = -1;
            I(ii+1) = 0;
        elseif (SOC(ii) < 0.01)
            holup = 1;
            I(ii+1) = 0;
        else
            I(ii+1) = I(ii);
        end
    end
    
    states = [c_c(:,ii+1);c_t(:,ii+1)];
   
    V(ii) = Vfun(states,[Q;I(ii)]) + 0*randn(1,1);
end

%%

figure(1)
subplot(2,1,1)
plot(c_c(1,:),'r')
hold on
plot(c_c(3,:),'--k')
legend('$c_{1a}$','$c_{2a}$')
subplot(2,1,2)
plot(c_c(2,:),'b')
hold on
plot(c0,'--r')
legend('$c_{1c}$','$c_{0c}$')
hold off
figure(2)
hold on
for ii = 1:size(c_t,1)
    plot(c_t(ii,:))
end
hold off
legend('$c_{1a}$','$c_{1c}$','$c_{2a}$')
figure(4)
subplot(3,1,1)
plot(V)
hold on
% plot(real(yp),'--r')
legend('Stack voltage','Estimated stack voltage')
subplot(3,1,2)
plot(I)
legend('Stack current')
subplot(3,1,3)
plot(SOC)
hold on
plot(realSOC,'--r')
hold off
legend('Anolyte-side SOC')

states = [c_c;c_t];
figure(5)
for ii = 1:size(states,1)
nexttile
plot(states(ii,:));
hold on
% plot(real(xp(ii,:)),'--r');
hold off
legend('Real concentration','Estimated concentration')
end