clc; clear; close all;

Data = xlsread('Ming-Li School 1999-9-21 data.xlsx');
data = Data(:,[8,15,21,28]);    % x: Time series
time_step = 0.005;
n = length(data);
fs = 1 / time_step;
frequencies = (0:n-1)*fs/n;
t = 0:0.005:102.395;

mean_data = zeros(n, 4);
for i = 1:4
    mean_data(:,i) = mean(data(:,i));
end
zeromean = data - mean_data;

%% Butterworth & FFT
forder = 10;
cff = 5;
[b,a] = butter(forder,cff/(fs/2),'low');
Ab = filtfilt(b,a,zeromean);
% figure
% freqz(b,a)
amplitude_spectrum = abs(fft(Ab))*2/n;
aa = [8,15,21,28];

figure
for i = 1:4
    subplot(4,1,i)
    plot(frequencies,amplitude_spectrum(:,i))
    xlim([0,10]); ylim([0,8]);
    xlabel('Frequency (Hz)','FontSize',11); ylabel('Amplitude','FontSize',11);
    title(['Floor ', num2str(aa(i)), ' Spectrum'],'FontSize',12);
    grid on
end
set(gcf,'unit','normalized','position',[0,0.2,0.42,0.5]) ;
saveas(gcf,fullfile(cd,['1_FFT spectrum.png']))

%% Spectrogram
F = 0:0.01:5.;
type_coef = 'energy';
sigma = 1.0;
NewY = zeromean;
for i = 1:4
    [ coefs ] = MCMW_VCF(NewY(:,i),F,fs,type_coef,sigma);
    AxisValue = double(max(max(abs(NewY))));
    figure
    subplot(3,3,[1,2])
    plot(t,NewY(:,i))
    ylabel('Accel. (gal)','FontSize',11)
    xlim([20 100]);ylim([-250 250]);
    title(['\bfAccel. Ch ', num2str(aa(i))],'FontSize',12)
    grid on
    
    subplot(3,3,[4,5,7,8])
    imagesc(t,F,coefs)
    set(gca,'YDir','normal');colormap Jet;
    axis([20 100 0 5])
    text(25,4.5,['\sigma = ',num2str(sigma)],'BackgroundColor',[1 1 1],'FontSize',12)
    xlabel('Time (sec)','FontSize',11);ylabel('Freq. (Hz)','FontSize',11);
    title(['\bfMCMW-VCF (Ch ', num2str(aa(i)),'),\sigma = ',num2str(sigma)],'FontSize',12)
    colorbar;
    grid on
    
    subplot(3,3,[6,9])
    MargSpec = sum(coefs,2);
    plot(MargSpec,F,'LineWidth',2)
    xlabel('Amplitude','FontSize',11)
    grid on

    set(gcf,'unit','normalized','position',[0,0.2,0.42,0.5]) ;
    saveas(gcf,fullfile(cd,['1_MCMW-VCF (Ch ', num2str(aa(i)),').png']))
end

% MCMW_VCF function
function [coefs] = MCMW_VCF(Acc,f,fs,type_coef,sigma)
    % Acc = input data
    % f = analysis frequency
    % fs = sampling rate
    % type_coef = output type of wavelet coefficients ('complex','absolute','energy')
    % sigma = the modified term of complex morlet wavelet
    % coefs= wavelet coefficients
    if size(Acc,2) == 1
        Acc = Acc';
    end
     Acc1=hilbert(Acc); 
     Acc=Acc1./abs(Acc1);
    scales = 1;
    i = (-1)^0.5;
    for s = 1:length(f)
         j = 1;
        for w = (0:length(Acc)-1)*fs/length(Acc)*2*pi
            psi(j) = scales^0.5*((2*pi)^0.5*sigma*exp(-0.5*(w*scales-2*pi*f(s))^2*sigma^2));
            j = j+1;
        end
        coefs(s,:) = 1/2/pi*(ifft(fft(Acc).*conj(psi)));
    end

    switch type_coef
        case 'complex'
            coefs = coefs;
        case 'absolute'
            coefs = abs(coefs);
        case 'energy'
            coefs = abs(coefs);
            coefs = coefs.*coefs;
    end
end