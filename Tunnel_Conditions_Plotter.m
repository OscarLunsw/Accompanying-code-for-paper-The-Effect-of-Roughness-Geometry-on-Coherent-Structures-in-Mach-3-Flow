clear all
close all
g = 1.4;
R = 287;
load('Cylinder_NearWake_Run31-Pt.mat')
pt = amplitude;
load('Cylinder_NearWake_Run31-Ps.mat')
ps = amplitude;
t1 = t;
load('Cylinder_NearWake_Run31-T')
t0 = t;
T = amplitude+273.15;

Ti = interp1(t0,T,t1);
M = sqrt(((pt./ps).^((g-1)/g)-1).*2/(g-1));

Ts = Ti.*(1+(g-1)/2.*M.^2).^-1;
Ue = M.*sqrt(g*R.*Ts);


subplot(4,1,1)
plot(t1,M,'k')
xlabel('Time, s')
ylabel('Mach Number')

subplot(4,1,2)
plot(t1,pt,'k')
xlabel('Time, s')
ylabel('P_0, Pa')


subplot(4,1,3)
plot(t1,Ti,'k')
xlabel('Time, s')
ylabel('T_0, K')


subplot(4,1,4)
plot(t1,Ue,'k')
xlabel('Time, s')
ylabel('U_\infty, m/s')

