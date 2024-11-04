clc; clear; close all;
data = xlsread('No.110103_X-dir (mean).xlsx');
t = data(:,1);
Base_acc = data(:,2)*0.01;
Top_acc = data(:,3)*0.01;
dt = t(2)-t(1);
SR = 1/dt;
%% (a)
base_acc = detrend(Base_acc);
top_acc = detrend(Top_acc);

Base_vel = cumtrapz(t,base_acc);
Top_vel = cumtrapz(t,top_acc);

base_vel = detrend(Base_vel);
top_vel = detrend(Top_vel);

Base_dis = cumtrapz(t,base_vel);
Top_dis = cumtrapz(t,top_vel);

top_dis = detrend(Top_dis,10);
base_dis = detrend(Base_dis,10);

rela_dis = top_dis - base_dis;
figure;
plot(t,top_dis);title('top displacement');xlabel('t');ylabel('m');
figure;
plot(t,base_dis);title('base displacement');xlabel('t');ylabel('m');
figure;
plot(t,rela_dis);title('(a) relative displacement');xlabel('t');ylabel('m');
%% (b)
Rela_acc = Top_acc - Base_acc;
rela_acc = detrend(Rela_acc);
Rela_vel = cumtrapz(t,rela_acc);
rela_vel = detrend(Rela_vel);
Rela_dis = cumtrapz(t,rela_vel);
rela_dis = detrend(rela_dis,10);
figure;
plot(t,rela_dis);title('(b) relative displacement');xlabel('t');ylabel('m');
%% (c)Ormsby
Rela_ACC = Top_acc - Base_acc;
fcut = 0.1;
ft = 0.15;
A = ormsby(Rela_ACC,ft,fcut,length(t),SR,dt);
A = A-mean(A(1:10*SR));
V = cumtrapz(t,A);
v = detrend(V);
d = cumtrapz(t,v);
figure
plot(t,d)
xlabel('t')
ylabel('m')
title('Relative Displacement(Ormsby)')
%% (c)Butterworth
forder = 10;
cff = 30;
[b,a] = butter(forder,cff/(SR/2),'low');
topacc = filtfilt(b,a,Top_acc);
baseacc = filtfilt(b,a,Base_acc);
relaacc = topacc-baseacc;
Relaacc = relaacc-mean(relaacc(1:10*SR));
Relavel = cumtrapz(t,Relaacc);
relavel = detrend(Relavel);
Reladis = cumtrapz(t,relavel);
figure;
plot(t,Reladis)
xlabel('t')
ylabel('m')
title('Relative Displacement(Butterworth)')
figure;
hold on
plot(t,d);
plot(t,Reladis);
xlabel('t');ylabel('m');title('Relative Displacement')
legend('Ormsby','Butterworth')
hold off


