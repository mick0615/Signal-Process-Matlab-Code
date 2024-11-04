clc;clear;close all;
data = xlsread('Ming-Li School 1999-9-21 data.xlsx');
base_acc = data(:,28);
roof_acc = data(:,21);
SR = 200 ;
fs = 200;
time = (0 : 1/SR : (length(data(:,1))-1)*1/SR)';
%% low pass and zero mean
%low pass 20HZ
order = 6;
cff = 20;
[b,a] = butter(order,cff/(SR/2),'low');
base_acc = filtfilt(b,a,base_acc);
roof_acc = filtfilt(b,a,roof_acc);
%zero mean
base_acc = base_acc-mean(base_acc);
roof_acc = roof_acc-mean(roof_acc);
%% plot wavelet coherence (remove arrow)
ch28 = base_acc;
ch21 = roof_acc;
[wcoh,~,f,coi] = wcoherence(ch28,ch21,fs);
figure;
surf(time,f,wcoh,'LineStyle','none','FaceLighting','phong');grid on
colormap('jet');
xlabel('Time (sec)');
ylabel('Freq. (Hz)');
xlim([0 100])
ylim([0 8])
view(0,90);colorbar;
title('wavelet coherence between basement response(ch28) and roof response(ch21)'...
    ,FontSize=15);