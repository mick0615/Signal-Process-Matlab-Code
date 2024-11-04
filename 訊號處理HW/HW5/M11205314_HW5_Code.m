%張原嘉-HW5-Code
clc; clear; close all;
data = importdata('HW-5-Pre-event Data of No-110103 NCREE.xlsx');
L   = length(data(:,1));
Fs  = 200; % (Hz)
T   = 1/Fs; % (s)
t   = (0:L-1)*T; % detla t
%% Plot Time Domain data for Floor 1-8
figure;
for i = 1:8
    subplot(4,2,i)
    plot(t, data(:,i))
    xlim([0 t(end)])
    ylim([-0.2 0.2])
    xlabel('Time (sec)')
    ylabel('Acceleration (g)')
    title( ['Floor ', num2str(i), ' in Time Domain' ] )
end
% Plot Time Domain data for Floor 9-14
figure;
for i = 1:6
    subplot(4,2,i)
    plot(t, data(:,i+8))
    xlim([0 t(end)])
    ylim([-0.2 0.2])
    xlabel('Time (sec)')
    ylabel('Acceleration (g)')
    title( ['Floor ', num2str(i+8), ' in Time Domain' ] )
end
%% Data Pre-Processing (Low Pass Filter)
cff = 20;
Acc = zeros(L,14);
[Para_B, Para_A] = butter(6,cff/(0.5*Fs), 'low');
for I = 1:size(data, 2)
    Acc(:, I) = filtfilt(Para_B, Para_A, data(:, I));
end
clear Para_A Para_B

% Zero Mean
for i = 1:14
    Acc(:, i) = Acc(:, i) - mean(data(:, i));
end
%% Auto-Spectral Density Function
cor_l = L/4; % correlation data length
f     = Fs*(0:(L-1)/2)'/L;

% create zero matrix after pre-processing
[~, ch] = size(Acc);
for i = 1: ch
    for j = 1:ch
        [Ryy(i,j,:)] = xcorr(Acc(:,i), Acc(:,j), cor_l, 'unbiased');   
    end
end

pt = size(Ryy, 3);
Ryy = Ryy(:, :, 1:round(pt/2));
pt = size(Ryy, 3);
Ryy = Ryy(: ,: , end: -1:1);
Ryy = Ryy(:, :, 1:round(pt/2));
pt = size(Ryy, 3);
Freq = (0: pt-1)/(pt/Fs);
PSD = zeros(ch, ch, pt);

for i = 1: ch
    for j = 1:ch
        Lu(i, j, :) = (Ryy(i, j, :) + Ryy(j, i, :)) /2;
        Qu(i, j, :) = (Ryy(i, j, :) - Ryy(j, i, :)) /2;   
        Lf(i, j, :) = real(fft(Lu(i, j, :)));
        Qf(i, j, :) = imag(fft(Qu(i, j, :)));
    end
end

PSD = Lf - (Qf*1i);
%% 2D PSD
for i = 1:1:size(PSD,3)
    PSD2D(:,i)=PSD(:,1,i);
end

PSD2D = PSD2D.';
% Plot 2D PSD
figure;
for i = 1:8
    subplot(4,2,i)
    plot(Freq, PSD2D(:,i))
    xlim([0 20])
    ylim([-1*10^-4 3*10^-4])
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
    title( ['Floor ', num2str(i), ' PSD' ] )
end

figure;
for i = 1:6
    subplot(4,2,i)
    plot(Freq, PSD2D(:,i+8))
    xlim([0 20])
    ylim([-1*10^-4 3*10^-4])
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
    title( ['Floor ', num2str(i+8), ' PSD' ] )
end
%% Single Value Decomposition
L_PSD = size(PSD, 3);

for i = 1: L_PSD
    [ U(:, :, i), SS, V(:, :, i) ] = svd(PSD(:, :, i).' ); % perform singular value decomposition
    Sdiag(:, i) = diag(SS);  % Extract Diagonal
    U1(:, i) = U(:, 1, i);% Extract Values
end
% plot SVD
figure;
semilogy(Freq, Sdiag(1, :), 'k', Freq, Sdiag(2, :)/5 , 'b',  Freq, Sdiag(3, :)/10, 'r')
xlim([0 20])
xlabel('Freqency (Hz)')
ylabel('Amplitude')
title( ['Singular Value Spectrum'] )
legend('1st SV', '2nd SV', '3rd SV')
%% Extract indexes for each principle frequency
for i = 1: size(PSD,3)
    if Freq(i) > 1.55944 && Freq(i) < 1.55946
        M1 = i;
    elseif Freq(i) > 3.50876 && Freq(i) < 3.50878
        M2 = i;
    elseif Freq(i) > 5.45808 && Freq(i) < 5.45810
        M3 = i;
    elseif Freq(i) > 7.40740 && Freq(i) < 7.40742
        M4 = i;
    elseif Freq(i) > 9.74658 && Freq(i) < 9.74660
        M5 = i;
    elseif Freq(i) > 13.2553 && Freq(i) < 13.2555
        M6 = i;
    elseif Freq(i) > 14.8147 && Freq(i) < 14.8149
        M7 = i;
    elseif Freq(i) > 17.5438 && Freq(i) < 17.5440
        M8 = i;
    elseif Freq(i) > 19.882 && Freq(i) < 19.884
        M9 = i;
    end
end
% Create empty mode shape array
fr = zeros(10, 9);
    
% Fill mode shape array at the identified indices
for i = 1:14
    fr(i+1,:) = [U1(i,M1),U1(i,M2),U1(i,M3),U1(i,M4),U1(i,M5),U1(i,M6),U1(i,M7),U1(i,M8),U1(i,M9)] ;
end
    
% create array of floor number
floorplot = (0:1:14)';

% creat vertical line to plot against
vertical = zeros(15,1);
 
% Radian and Angle
R = abs(fr);
theta = angle(fr);

% %% Plot Mode Shape and Pole Location
% for j = 1:9
%     figure()
%     % mode shape
%     subplot(1,2,1)
%     plot(vertical, floorplot, '--r', fr(:, j), floorplot, '-bo')
%     xlim([-1 1])
%     ylabel('Floor Number')
%     xlabel('Amplitude')
%     title( ['Mode Shape ', num2str(j)])
%     legend('Floor plot', 'Mode shape')
%     
%     % location of poles
%     subplot(1,2,2)
%     polarplot(theta(:, j), R(:, j), 'ro');
%     rlim([0 1])
%     title('Location of poles')
% end