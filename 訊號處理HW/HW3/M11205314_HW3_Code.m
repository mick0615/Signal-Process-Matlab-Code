%張原嘉-HW3-Code
clear, clc, close all;
data = xlsread('No.110103_X-dir (mean).xlsx');
t = data(:,1);
base_x_ch1 = data(:,2);
roof_x_ch19 = data(:,3);
L = length(base_x_ch1);
del_t = t(2);
fs=1/del_t;
L_time = (L-1) * del_t;
% Plot acceleration data of base and roof in time domain
figure
subplot(211),plot(t,base_x_ch1)
xlim([0,L_time])
xlabel('Time (s)')
ylabel('Acceleration (gal)')
title('Acceleration vs Time for Basement, Global X Direction')
subplot(212),plot(t,roof_x_ch19)
xlim([0,L_time])
xlabel('Time (s)')
ylabel('Acceleration (gal)')
title('Acceleration vs Time for 14th Floor, Global X Direction')
% Calculate Auto-Correlation of Basement Data in mothod 2
[auto_base_m2,lag_m2] = xcorr(base_x_ch1,length(base_x_ch1)/2,'unbiased');
start = find(lag_m2 == 0);
autobase = auto_base_m2(start:end);
lagbase = lag_m2(start:end);
% Plot the 2 steps of Auto-Correlation Basement Data in mothod 2
figure
sgtitle('Method 2: Auto-Correlation Calculation')
subplot(211),plot(lag_m2,auto_base_m2)
xlim([lag_m2(1),lag_m2(end)])
xlabel('Point')
ylabel('Amplitude')
title('Auto-Correlation Basement, L/2 Overlap (Rxx)')
grid on
subplot(212),plot(lagbase,autobase)
xlim([lag_m2(1),lag_m2(end)])
xlabel('Point')
ylabel('Amplitude')
title('Auto-Correlation Basement,Second Half of Data (Rxx)')
grid on
% Calculate Auto-Correlation of Roof Data

[auto_roof,lag_roof] = xcorr(roof_x_ch19,length(roof_x_ch19)/2,'unbiased');
start = find(lag_roof == 0);
autoroof = auto_roof(start:end);
lagroof = lag_roof(start:end);
% Plot for Auto-Correlation of Roof Data
figure
plot(lagroof,autoroof,'LineWidth',1.5)
xlim([lagroof(1),lagroof(end)])
xlabel('Point')
ylabel('Amplitude')
title('Auto-Correlation of Roof Acceleration Data (Ryy)')
grid on
% Calculate Cross-Correlation of Basement and Roof Data
[cross_xy,lag_xy] = xcorr(base_x_ch1,roof_x_ch19,round(length(base_x_ch1)/2),'unbiased');
cross_xy_flip = flip(cross_xy(1:find(lag_xy==0)));
lag_xy_2 = lag_xy(find(lag_xy==0):end);
%Plot results for Rxy cross-correlation
figure
subplot(211)
plot(lag_xy_2,cross_xy_flip,'LineWidth',1.5)
xlim([lag_xy_2(1),lag_xy_2(end)])
xlabel('Point')
ylabel('Amplitude')
title('Cross-Correlation of Basement and Roof Acceleration Data (Rxy)')
grid on
[cross_yx,lag_yx] = xcorr(roof_x_ch19,base_x_ch1,round(length(base_x_ch1)/2),'unbiased');
cross_yx_flip = flip(cross_yx(1:find(lag_yx==0)));
lag_yx_2 = lag_yx(find(lag_yx==0):end);
%Plot results for Ryx cross-correlation
subplot(212)
plot(lag_yx_2,cross_yx_flip,'LineWidth',1.5)
xlim([lag_yx_2(1),lag_yx_2(end)])
xlabel('Point')
ylabel('Amplitude')
title('Cross-Correlation of Basement and Roof Acceleration Data (Ryx)')
grid on
%Find time lag at the maximum amplitude of the cross-correlation function
peak_Rxy = max(cross_xy_flip);
peak_Ryx = max(cross_yx_flip);
peak_lag_point = lag_xy_2(find(cross_xy_flip==peak_Rxy));
lag_time_x2y = del_t*peak_lag_point
height = 50;
shear_wave_propagation = height / lag_time_x2y
disp(['shear_wave_propagation = ',num2str(shear_wave_propagation),'m/s'])