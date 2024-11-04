%張原嘉 M11205314 HW9 Code
clc;clear;close all;
data = xlsread('HW9-data.xlsx');
F = 0:0.1:8;
Fs = 200;
type_coef = 'energy';
RCF4_A4 = data(:,7);
RCF4_A10 = data(:,6);
F4_acc = RCF4_A4 - RCF4_A10;
t = data(:,1);

%low pass 20Hz Freuency 
order = 6;
cff = 20;
[b,a] = butter(order,cff/(Fs/2),'low');
F4_acc = filtfilt(b,a,F4_acc);
% down sampling to 50Hz
SR = 1/200;
NewSR = 1/50;
F4_acc = F4_acc([1:NewSR/SR:end]);
t = t([1:NewSR/SR:end]);
New_acc = F4_acc-mean(F4_acc);

N = length(t);
% nfft = 2^nextpow2(N);
% Df=1/(nfft*SR);
faxis =linspace(0,128,N/2+1);
faxis= transpose(faxis);
Signal_FFT = abs(fft(New_acc , N));
plot(faxis,Signal_FFT(1:N/2+1)) ; grid on
axis([ 0 10 0 max(Signal_FFT(1:N/2+1))*1.1])
title('FFT','FontSize',12)
ylabel('Amplitude') ; xlabel('Freq.(Hz)')
ylabel('Amplitude','FontSize',11) ; ylabel('Freq.(Hz)','FontSize',11)
set(gcf,'unit','normalized','position',[0.3,0.3,0.4,0.5]);

dt = 1/50;
fs = 50; % sampling rate
Ndata = length(New_acc); % number of points
%% Sigma = 2
sigma = 2;
[ coefs ] = MCMW_VCF( New_acc , F , fs , type_coef , sigma);
AxisValue = double(max(max(abs(New_acc))));

figure (4)
subplot(3,3,[1,2])
plot(t , New_acc) ;  ylabel('\bfAccel. (g)') ; grid on
axis([0 40 -AxisValue*1.1 AxisValue*1.1]) ;  %set(gca,'fontsize',10)
xlim([10 90]) ; xticks(0:5:90);
set(gca , 'XTickLabelRotation',0)
title(['\bfnewData'],'FontSize',12)

subplot(3,3,[4,5,7,8])
imagesc(t,F,coefs) ; grid on ; set(gca,'YDir','normal') ; colormap Jet ;
axis([10 90 0 8]) ; set(gca , 'XTickLabelRotation',0) ; xticks(0:10:90);
text(15, 7.5,['\sigma = ',num2str(sigma)],'BackgroundColor',[1 1 1],'FontSize',12)
xlabel('\bfTime (sec)') ; ylabel('\bfFreq. (Hz)') ;
colorbar('position',[0.63 0.11 0.023 0.515],'FontSize',8)

% Marginal Spectrum
subplot(3,3,[6,9])
MargSpec = sum(coefs,2);
plot(MargSpec , F ,'LineWidth',1);
xlabel('Amplitude','FontSize',11)

set(gcf,'unit','normalized','position',[0.25,0.25,0.42,0.5]) ;
sgtitle('Sigma = 2') ;
%[ coefs ] = MCMW_VCF( NewY , F , fs , type_coef , sigma);
%% sigma = 1
sigma = 1;
[ coefs ] = MCMW_VCF( New_acc , F , fs , type_coef , sigma);
AxisValue = double(max(max(abs(New_acc))));

figure;
subplot(3,3,[1,2])
plot(t , New_acc) ;  ylabel('\bfAccel. (g)') ; grid on
axis([0 40 -AxisValue*1.1 AxisValue*1.1]) ;  %set(gca,'fontsize',10)
xlim([10 90]) ; xticks(0:5:90);
set(gca , 'XTickLabelRotation',0)
title(['\bfnewData'],'FontSize',12)

subplot(3,3,[4,5,7,8])
imagesc(t,F,coefs) ; grid on ; set(gca,'YDir','normal') ; colormap Jet ;
axis([10 90 0 8]) ; set(gca , 'XTickLabelRotation',0) ; xticks(0:10:90);
text(15, 7.5,['\sigma = ',num2str(sigma)],'BackgroundColor',[1 1 1],'FontSize',12)
xlabel('\bfTime (sec)') ; ylabel('\bfFreq. (Hz)') ;
colorbar('position',[0.63 0.11 0.023 0.515],'FontSize',8)

% Marginal Spectrum
subplot(3,3,[6,9])
MargSpec = sum(coefs,2);
plot(MargSpec , F ,'LineWidth',1);
xlabel('Amplitude','FontSize',11)

set(gcf,'unit','normalized','position',[0.25,0.25,0.42,0.5]) ;
sgtitle('Sigma = 1') ;
%% sigma = 0.2
sigma = 0.2;
[ coefs ] = MCMW_VCF( New_acc , F , fs , type_coef , sigma);
AxisValue = double(max(max(abs(New_acc))));

figure;
subplot(3,3,[1,2])
plot(t , New_acc) ;  ylabel('\bfAccel. (g)') ; grid on
axis([0 40 -AxisValue*1.1 AxisValue*1.1]) ;  %set(gca,'fontsize',10)
xlim([10 90]) ; xticks(0:5:90);
set(gca , 'XTickLabelRotation',0)
title(['\bfnewData'],'FontSize',12)

subplot(3,3,[4,5,7,8])
imagesc(t,F,coefs) ; grid on ; set(gca,'YDir','normal') ; colormap Jet ;
axis([10 90 0 8]) ; set(gca , 'XTickLabelRotation',0) ; xticks(0:10:90);
text(15, 7.5,['\sigma = ',num2str(sigma)],'BackgroundColor',[1 1 1],'FontSize',12)
xlabel('\bfTime (sec)') ; ylabel('\bfFreq. (Hz)') ;
colorbar('position',[0.63 0.11 0.023 0.515],'FontSize',8)

%Marginal Spectrum
subplot(3,3,[6,9])
MargSpec = sum(coefs,2);
plot(MargSpec , F ,'LineWidth',1);
xlabel('Amplitude','FontSize',11)

set(gcf,'unit','normalized','position',[0.25,0.25,0.42,0.5]) ;
sgtitle('Sigma = 0.2') ;
%% 9-2
data2 = xlsread('3Story.xlsx');
A19 = data2(:,2);
A24 = data2(:,3);
%low pass 20Hz Freuency and zero mean
order = 6;
cff = 20;
[b,a] = butter(order,cff/(Fs/2),'low');
t = data2(:,1);
A19 = filtfilt(b,a,A19);
A24 = filtfilt(b,a,A24);
A19 = A19-mean(A19);
A24 = A24-mean(A24);
[wcoh,~,f,coi] = wcoherence(A19,A24,Fs);
figure;
surf(t,f,wcoh,'LineStyle','none','FaceLighting','phong');grid on
colormap('jet');
xlabel('Time (sec)');
ylabel('Freq. (Hz)');
ylim([0 70])
view(0,90);colorbar;
title('wavelet coherence between A19 and A24',FontSize=15);

function [ coefs ] = MCMW_VCF( Acc,f,fs,type_coef,sigma)
% Acc = input data
% f = analysis frequency
% fs = sampling rate
% 
% type_coef = output type of wavelet coefficients ('complex', 'absolute', 'energy')
% sigma = the modified term of complex morlet wavelet
% coefs= wavelet coefficients

% -----------------------------------------------------------------------

if size(Acc,2) == 1
    Acc = Acc';
end
% Acc1=hilbert(Acc); 
% Acc=Acc1./abs(Acc1);

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


