%張原嘉-HW4-Code
clear all;clc
data = xlsread('No.110103_X-dir (mean).xlsx');
t = data(:,1);
x = data(:,2);
y = data(:,3);
ts = t(2);
fs = 1/ts;
cor_l = 1024;
% input(x) auto spectrum density function(via correlation function)
[Rxx,lag_input] = xcorr(x,cor_l,'unbiased');
pt = round((size(Rxx,1)/2));
Rxx = Rxx(1:pt);
Rxx = Rxx(end:-1:1);
pt = round((size(Rxx,1)/2));
Rxx = Rxx(1:pt);
Sxx = abs(fft(Rxx));
fxx = (0:pt-1)*(((1/fs)*pt)^-1);
figure;
subplot(2,1,1);
plot(fxx(1:pt/2),Sxx(1:pt/2)*2/pt);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('via Correlation Function','\bf Auto-spectrum Density Function (Sxx)');
xlim([0,5])
% input(x) auto spectrum density function(via finite fourier transform)
subplot(2,1,2)
pt2=round(size(x,1));
T=size(x,1)*ts;
X=fft(x);
Xi=conj(fft(x));
SXX=(X.*Xi)/T;
FXX = (0:pt2-1)*(((1/fs)*pt2)^-1);
plot(FXX,SXX);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('via direct fourier transform','\bf Auto-spectrum Density Function (Sxx)');
xlim([0,5]);

% output(y) auto spectrum density function(via correlation function)
[Ryy,lag_output] = xcorr(y,cor_l,'unbiased');
pt = round((size(Ryy,1)/2));
Ryy = Ryy(1:pt);
Ryy = Ryy(end:-1:1);
pt = round((size(Ryy,1)/2));
Ryy = Ryy(1:pt);
Syy = abs(fft(Ryy));
fyy = (0:pt-1)*(((1/fs)*pt)^-1);
figure;
subplot(2,1,1);
plot(fyy(1:pt/2),Syy(1:pt/2)*2/pt);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('via Correlation Function','\bf Auto-spectrum Density Function (Syy)');
xlim([0,5])
% output(y) auto spectrum density function(via finite fourier transform)
subplot(2,1,2);
pt2=round(size(y,1));
T=size(y,1)*ts;
Y=fft(y);
Yi=conj(fft(y));
SYY=(Y.*Yi)/T;
FYY = (0:pt2-1)*(((1/fs)*pt2)^-1);
plot(FYY,SYY);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('via direct fourier transform','\bf Auto-spectrum Density Function (Syy)');
xlim([0,5]);

% input(x) and output(y) cross-spectrum density function(via correlation)
[Rxy,lag_xy] = xcorr(x,y,cor_l,'unbiased');
[Ryx,lag_yx] = xcorr(y,x,'unbiased');
pt = round((size(Rxy,1)/2));
Rxy = Rxy(1:pt);
Ryx = Ryx(1:pt);

Rxy = Rxy(end:-1:1);
Ryx = Ryx(end:-1:1);

pt = round((size(Rxy,1)/2));
Rxy = Rxy(1:pt);
Ryx = Ryx(1:pt);
lxy = 0.5*(Rxy+Ryx);
qxy = 0.5*(Rxy-Ryx);
Lxy = real(fft(lxy));
Qxy = imag(fft(qxy));
Sxy = Lxy-Qxy*1i;
fxy = (0:pt-1)*(((1/fs)*pt)^-1);
figure;
subplot(2,1,1);
plot(fxy(1:pt/2),Sxy(1:pt/2)*2/pt);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('via correlation function','\bf Auto-spectrum Density Function (Sxy)');
xlim([0,5])
% input(x) and output(y) cross-spectrum density function(via direct fourier)
subplot(2,1,2);
pt2=round(size(x,1));
T=size(x,1)*ts;
X=fft(x);Y=fft(y);
SXY=(Xi.*Y)/T;
FXY = (0:pt2-1)*(((1/fs)*pt2)^-1);
plot(FXY,SXY);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('via direct fourier transform','\bf Auto-spectrum Density Function (Sxy)');
xlim([0,5]);
%Calculate the coherence function between input(x) and output(x)
figure;
subplot(2,1,1);
gamma1 = (abs(Sxy)).^2/(Sxx.*Syy);
plot(fxy,gamma1);
xlim([0 5]);
xlabel('Frequency(Hz)');
title('Coherence (via Correlation Function)');
subplot(2,1,2);
gamma2 = (abs(SXY)).^2/(SXX.*SYY);
plot(FXY,gamma2);
xlim([0 5]);
xlabel('Frequency(Hz)');
title('Coherence (via direct fourier transform)');