clc;clear;close all;
data = xlsread('Data set_1story_2bay_RCF-Four_ Specimen_Test_Data.xlsx');
RCF6_A4 = data(:,3);
RCF6_A10 = data(:,2);
t = data(:,1);
N = length(t);
SR = 200; % sampling rate, Hz
T = 1/SR; % sampling period, sec
dt = 1/SR;
u = RCF6_A4;
Level = 9;
wname = 'bior6.8';
index = 0;
for jlevel=1:Level
    index = [index fliplr(index)+2^(jlevel-1)];
end
index = index+1;    
% Plot wavelet packet tree wpt
WPTree = wpdec(u,Level,wname);
%plot(WPTree)

Amp = zeros(length(u),2^Level);
for icomp = 1 : 1 : 2^Level
     RconstCoef(:,icomp)= wprcoef(WPTree,[Level,index(icomp)-1]);
     Amp(:,icomp) = abs(hilbert(RconstCoef(:,icomp)));
end
Energy_bior_9 = Amp.*Amp;
Fvec = linspace(0,1/dt/2,2^Level)';
Fint = find(Fvec<50);
% Plot time-freq spectrogram
figure('Position',[200 50 500 300])
surf(t,Fvec(Fint),Energy_bior_9(:,Fint)','LineStyle','none','FaceColor','interp','FaceLighting','phong','EdgeColor','interp')
view(0,90),colormap('jet'),ylim([0 8]),xlabel('Time (sec)'),ylabel('Frequency (Hz)'),xlim([0 t(end)]);
colorbar
title("Original WPT,Level = 9 (bior6.8)")
EnergySTORE_bior_9(:,:,1) = Energy_bior_9;
%% Enhanced TF representation(select first 1 singular value to reconstruct the spectrogram)
[U S V] = svd(Energy_bior_9);
k = 1;
New_S = S(1:k,1:k);
New_Energy_9 = U(:,1:k) * New_S * V(:,1:k)';
% Plot time-freq spectrogram
figure('Position',[200 50 500 300])
surf(t,Fvec(Fint),New_Energy_9(:,Fint)','LineStyle','none','FaceColor','interp','FaceLighting','phong','EdgeColor','interp')
view(0,90),colormap('jet'),ylim([0 8]),xlabel('Time (sec)'),ylabel('Frequency (Hz)'),xlim([0 t(end)]);
colorbar
title("Enhanced TF WPT,Level = 9 (bior6.8)")