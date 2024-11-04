%張原嘉-HW6-Code
clc;clear;close all;
data = importdata('HW-5-Pre-event Data of No-110103 NCREE.xlsx');
SR = 200;
Fs = SR;
dt = 1/SR;
%L = length(data(:,1));
L = 4096;
N = L;
t = ((0:L-1)*dt)';
for i = 1:7
    A(:,i) = data(:,i*2);
end
%% plot 2 4 6 8 floor(Time domain)
figure;
for i = 1:4
    subplot(4,1,i)
    plot(t,A(:,i));xlabel('time(sec)');ylabel('Acc(m/s2)');xlim([0 t(end)]);
    title(['Floor ',num2str(i*2),' Acceleration'])
end
%% plot 10 12 RF(Time domain)
figure;
for i = 1:3
    subplot(3,1,i)
    plot(t,A(:,4+i));xlabel('time(sec)');ylabel('Acc(m/s2)');xlim([0 t(end)]);
    title(['Floor ',num2str((i+4)*2),' Acceleration'])
end
%% low pass and zero mean
cff = 20;
[para_a para_b] = butter(6,cff/(0.5*SR),'low');
for i = 1:size(A,2)
    A(:,i) = filtfilt(para_a,para_b,A(:,i));
end
%zero pass
for i = 1:size(A,2)
    A(:,i) = A(:,i) - mean(A(:,i));
end
%% plot fourier amplitude 2 4 6 8 floor
figure;
for i = 1:4
    ffA = fft(A(:,i));
    mag = abs(ffA);
    fi = (0:L-1)*SR/L;
    subplot(4,1,i)
    plot(fi(1:L/2),mag(1:N/2)*2/N);xlabel('Hz');ylabel('AMP')
    title(['Floor ',num2str(i*2),' Frequency damain'])
end
%% plot fourier amplitude 10 12 RF
figure;
for i = 1:3
    ffA = fft(A(:,i+4));
    mag = abs(ffA);
    fi = (0:L-1)*SR/L;
    subplot(3,1,i)
    plot(fi(1:L/2),mag(1:L/2)*2/N);xlabel('Hz');ylabel('AMP')
    title(['Floor ',num2str((i+4)*2),' in Frequency domain'])
end
%% plot the distribution of the singular value
[pt,ch]=size(A);
totalpair = 7;
nrh = 500;
nch = pt-nrh+1;
X = zeros(ch*nrh,nch);
for i=1:ch
    X(i:ch:end,:)=hankel(A(1:nrh,i),A(nrh:pt,i));
end

if true
    [U,S,V]=svd(X,'econ');
else
    Scov=X*X';
    [U,S,V]=svd(Scov,'econ');
    clear Scov
end
S = diag(S);
figure;
plot(1:length(S),S,'-o','LineWidth',1.5);xlim([0,20]);
xlabel('Singular Value');ylabel('Singular Value');title('Singular Value Plot')
figure;
plot(1:length(S),100*cumsum(S)/sum(S),'-o','LineWidth',1.5);xlim([0 20]);
xlabel('Number of Singular Value');ylabel('Percentage of Accumulated Singular Values');
title('Singular Value Distribution,L=100')
%% select pairs of SV to reconstruct data
totalpair = 7;
X_recon = zeros([size(X),totalpair]);
for i=1:totalpair
    tmpidx = [2*i-1,2*i];
    if true
        tmp = U(:,tmpidx)*diag(S(tmpidx))*V(:,tmpidx)';
    else
        tempv1 = X'*U(:,tmpidx(1))/sqrt(S(tmpidx(1)));
        tempv2 = X'*U(:,tmpidx(2))/sqrt(S(tmpidx(2)));
        tmp = U(:,tmpidx)*diag(sqrt(S(tmpidx)))*[tempv1,tempv2];
    end
    X_recon(:,:,i)=tmp;
end
clear tmpidx tmp tmpv1 tmpv2
%% Reconstruction data
x_recon = zeros(pt,ch,totalpair);
for k=1:totalpair
    for i=1:ch
        tmph=X_recon(i:ch:end,:,k);
        tmph=flipud(tmph);
        for j=1:pt
            tmp=diag(tmph,j-nrh);
            x_recon(j,i,k)=mean(tmp);
        end
    end
end
clear tmph
ReCon_1=x_recon(:,1,1)+x_recon(:,1,2)+x_recon(:,1,3)+x_recon(:,1,4)+x_recon(:,1,5)+x_recon(:,1,6)+x_recon(:,1,7);
ReCon_2=x_recon(:,2,1)+x_recon(:,2,2)+x_recon(:,2,3)+x_recon(:,2,4)+x_recon(:,2,5)+x_recon(:,2,6)+x_recon(:,2,7);
ReCon_3=x_recon(:,3,1)+x_recon(:,3,2)+x_recon(:,3,3)+x_recon(:,3,4)+x_recon(:,3,5)+x_recon(:,3,6)+x_recon(:,3,7);
ReCon_4=x_recon(:,4,1)+x_recon(:,4,2)+x_recon(:,4,3)+x_recon(:,4,4)+x_recon(:,4,5)+x_recon(:,4,6)+x_recon(:,4,7);
ReCon_5=x_recon(:,5,1)+x_recon(:,5,2)+x_recon(:,5,3)+x_recon(:,5,4)+x_recon(:,5,5)+x_recon(:,5,6)+x_recon(:,5,7);
ReCon_6=x_recon(:,6,1)+x_recon(:,6,2)+x_recon(:,6,3)+x_recon(:,6,4)+x_recon(:,6,5)+x_recon(:,6,6)+x_recon(:,6,7);
ReCon_7=x_recon(:,7,1)+x_recon(:,7,2)+x_recon(:,7,3)+x_recon(:,7,4)+x_recon(:,7,5)+x_recon(:,7,6)+x_recon(:,7,7);
%% Floor 2 Response
figure;
subplot(8,1,1)
plot(t,A(:,1),'LineWidth',1.5);xlim([0 20])
hold on
plot(t,ReCon_1,'LineWidth',1.5);xlim([0 20])
xlabel('Time(sec)');ylabel('Resp.');ylim([-0.01 0.01]);legend('Original','Reconstructed')
for i = 1:totalpair
    subplot(8,1,i+1)
    plot(t,x_recon(:,1,i),'LineWidth',1.5);xlim([0 20])
    xlabel('Time(s)');ylabel('Resp.');ylim([-0.01 0.01]);legend(['Group',num2str(i)])
end
sgtitle('Floor 2 Response')
%% Floor 4 Response
figure;
subplot(8,1,1)
plot(t,A(:,2),'LineWidth',1.5);xlim([0 20])
hold on
plot(t,ReCon_2,'LineWidth',1.5);xlim([0 20])
xlabel('Time(sec)');ylabel('Resp.');ylim([-0.01 0.01]);legend('Original','Reconstructed')
for i = 1:totalpair
    subplot(8,1,i+1)
    plot(t,x_recon(:,2,i),'LineWidth',1.5);xlim([0 20])
    xlabel('Time(sec)');ylabel('Resp.');ylim([-0.01 0.01]);legend(['Group',num2str(i)])
end
sgtitle('Floor 4 Response')
%% Floor 6 Response
figure;
subplot(8,1,1)
plot(t,A(:,3),'LineWidth',1.5);xlim([0 20])
hold on
plot(t,ReCon_3,'LineWidth',1.5);xlim([0 20])
xlabel('Time(sec)');ylabel('Resp.');ylim([-0.01 0.01]);legend('Original','Reconstructed')
for i = 1:totalpair
    subplot(8,1,i+1)
    plot(t,x_recon(:,3,i),'LineWidth',1.5);xlim([0 20])
    xlabel('Time(sec)');ylabel('Resp.');ylim([-0.01 0.01]);legend(['Group',num2str(i)])
end
sgtitle('Floor 6 Response')
%% Floor 8 Response
figure;
subplot(8,1,1)
plot(t,A(:,4),'LineWidth',1.5);xlim([0 20])
hold on
plot(t,ReCon_4,'LineWidth',1.5);xlim([0 20])
xlabel('Time(sec)');ylabel('Resp.');ylim([-0.01 0.01]);legend('Original','Reconstructed')
for i = 1:totalpair
    subplot(8,1,i+1)
    plot(t,x_recon(:,4,i),'LineWidth',1.5);xlim([0 20])
    xlabel('Time(sec)');ylabel('Resp.');ylim([-0.01 0.01]);legend(['Group',num2str(i)])
end
sgtitle('Floor 8 Response')
%% Floor 10 Response
figure;
subplot(8,1,1)
plot(t,A(:,5),'LineWidth',1.5);xlim([0 20])
hold on
plot(t,ReCon_5,'LineWidth',1.5);xlim([0 20])
xlabel('Time(sec)');ylabel('Resp.');ylim([-0.01 0.01]);legend('Original','Reconstructed')
for i = 1:totalpair
    subplot(8,1,i+1)
    plot(t,x_recon(:,5,i),'LineWidth',1.5);xlim([0 20])
    xlabel('Time(sec)');ylabel('Resp.');ylim([-0.01 0.01]);legend(['Group',num2str(i)])
end
sgtitle('Floor 5 Response')
%% Floor 12 Response
figure;
subplot(8,1,1)
plot(t,A(:,6),'LineWidth',1.5);xlim([0 20])
hold on
plot(t,ReCon_6,'LineWidth',1.5);xlim([0 20])
xlabel('Time(sec)');ylabel('Resp.');ylim([-0.01 0.01]);legend('Original','Reconstructed')
for i = 1:totalpair
    subplot(8,1,i+1)
    plot(t,x_recon(:,6,i),'LineWidth',1.5);xlim([0 20])
    xlabel('Time(sec)');ylabel('Resp.');ylim([-0.01 0.01]);legend(['Group',num2str(i)])
end
sgtitle('Floor 12 Response')
%% RF Response
figure;
subplot(8,1,1)
plot(t,A(:,7),'LineWidth',1.5);xlim([0 20])
hold on
plot(t,ReCon_7,'LineWidth',1.5);xlim([0 20])
xlabel('Time(sec)');ylabel('Resp.');ylim([-0.01 0.01]);legend('Original','Reconstructed')
for i = 1:totalpair
    subplot(8,1,i+1)
    plot(t,x_recon(:,7,i),'LineWidth',1.5);xlim([0 20])
    xlabel('Time(sec)');ylabel('Resp.');ylim([-0.01 0.01]);legend(['Group',num2str(i)])
end
sgtitle('RF Response')
%% Plot FRF Floor 2
f = t-t(1);
f = Fs*f/f(end);
figure;
subplot(ch+1,1,1)
plot(f,2*abs(fft(A(:,1))),'--','LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_1)),'LineWidth',1.5)
xlabel('Frequency(Hz)');ylabel('Mag.');xlim([0 100]);legend('Original','Reconstructed')
for i=1:totalpair
    subplot(ch+1,1,i+1)
    plot(f,2*abs(fft(A(:,1))),'--','LineWidth',1.5)
    hold on
    plot(f,2*abs(fft(x_recon(:,1,i))),'LineWidth',1.5)
    xlabel('Frequency(Hz)');ylabel('Mag.');xlim([0 100]);
    legend('Original',['Group',sprintf('%d',i)])
end
sgtitle('Floor 2:FFT of SV components')
%% Plot FRF Floor 4
f = t-t(1);
f = Fs*f/f(end);
figure;
subplot(ch+1,1,1)
plot(f,2*abs(fft(A(:,2))),'--','LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_2)),'LineWidth',1.5)
xlabel('Frequency(Hz)');ylabel('Mag.');xlim([0 100]);legend('Original','Reconstructed')
for i=1:totalpair
    subplot(ch+1,1,i+1)
    plot(f,2*abs(fft(A(:,2))),'--','LineWidth',1.5)
    hold on
    plot(f,2*abs(fft(x_recon(:,2,i))),'LineWidth',1.5)
    xlabel('Frequency(Hz)');ylabel('Mag.');xlim([0 100]);
    legend('Original',['Group',sprintf('%d',i)])
end
sgtitle('Floor 4:FFT of SV components')
%% Plot FRF Floor 6
f = t-t(1);
f = Fs*f/f(end);
figure;
subplot(ch+1,1,1)
plot(f,2*abs(fft(A(:,3))),'--','LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_3)),'LineWidth',1.5)
xlabel('Frequency(Hz)');ylabel('Mag.');xlim([0 100]);legend('Original','Reconstructed')
for i=1:totalpair
    subplot(ch+1,1,i+1)
    plot(f,2*abs(fft(A(:,3))),'--','LineWidth',1.5)
    hold on
    plot(f,2*abs(fft(x_recon(:,3,i))),'LineWidth',1.5)
    xlabel('Frequency(Hz)');ylabel('Mag.');xlim([0 100]);
    legend('Original',['Group',sprintf('%d',i)])
end
sgtitle('Floor 6:FFT of SV components')
%% Plot FRF Floor 8
f = t-t(1);
f = Fs*f/f(end);
figure;
subplot(ch+1,1,1)
plot(f,2*abs(fft(A(:,4))),'--','LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_4)),'LineWidth',1.5)
xlabel('Frequency(Hz)');ylabel('Mag.');xlim([0 100]);legend('Original','Reconstructed')
for i=1:totalpair
    subplot(ch+1,1,i+1)
    plot(f,2*abs(fft(A(:,4))),'--','LineWidth',1.5)
    hold on
    plot(f,2*abs(fft(x_recon(:,4,i))),'LineWidth',1.5)
    xlabel('Frequency(Hz)');ylabel('Mag.');xlim([0 100]);
    legend('Original',['Group',sprintf('%d',i)])
end
sgtitle('Floor 8:FFT of SV components')
%% Plot FRF Floor 10
f = t-t(1);
f = Fs*f/f(end);
figure;
subplot(ch+1,1,1)
plot(f,2*abs(fft(A(:,5))),'--','LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_5)),'LineWidth',1.5)
xlabel('Frequency(Hz)');ylabel('Mag.');xlim([0 100]);legend('Original','Reconstructed')
for i=1:totalpair
    subplot(ch+1,1,i+1)
    plot(f,2*abs(fft(A(:,5))),'--','LineWidth',1.5)
    hold on
    plot(f,2*abs(fft(x_recon(:,5,i))),'LineWidth',1.5)
    xlabel('Frequency(Hz)');ylabel('Mag.');xlim([0 100]);
    legend('Original',['Group',sprintf('%d',i)])
end
sgtitle('Floor 10:FFT of SV components')
%% Plot FRF Floor 12
f = t-t(1);
f = Fs*f/f(end);
figure;
subplot(ch+1,1,1)
plot(f,2*abs(fft(A(:,6))),'--','LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_6)),'LineWidth',1.5)
xlabel('Frequency(Hz)');ylabel('Mag.');xlim([0 100]);legend('Original','Reconstructed')
for i=1:totalpair
    subplot(ch+1,1,i+1)
    plot(f,2*abs(fft(A(:,6))),'--','LineWidth',1.5)
    hold on
    plot(f,2*abs(fft(x_recon(:,6,i))),'LineWidth',1.5)
    xlabel('Frequency(Hz)');ylabel('Mag.');xlim([0 100]);
    legend('Original',['Group',sprintf('%d',i)])
end
sgtitle('Floor 12:FFT of SV components')
%% Plot FRF Floor RF
f = t-t(1);
f = Fs*f/f(end);
figure;
subplot(ch+1,1,1)
plot(f,2*abs(fft(A(:,7))),'--','LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_7)),'LineWidth',1.5)
xlabel('Frequency(Hz)');ylabel('Mag.');xlim([0 100]);legend('Original','Reconstructed')
for i=1:totalpair
    subplot(ch+1,1,i+1)
    plot(f,2*abs(fft(A(:,7))),'--','LineWidth',1.5)
    hold on
    plot(f,2*abs(fft(x_recon(:,7,i))),'LineWidth',1.5)
    xlabel('Frequency(Hz)');ylabel('Mag.');xlim([0 100]);
    legend('Original',['Group',sprintf('%d',i)])
end
sgtitle('Floor RF:FFT of SV components')
%% Plot FRF Floor RF
figure;
f = (0:L-1)*SR/L;
mag = abs(fft(A(:,7)));
subplot(ch+1,1,1)
hold on
plot(f(1:L/2),mag(1:L/2)*2/L,'b--','LineWidth',1.5)
maga = abs(fft(ReCon_7));
plot(f(1:L/2),maga(1:L/2)*2/L,'r','LineWidth',1.5)
xlabel('Frequency(Hz)');ylabel('Mag.');xlim([0 100]);legend('Original','Reconstructed')
hold off
for i = 1:totalpair
    subplot(ch+1,1,i+1)
    hold on
    mag = abs(fft(A(:,7)));
    magx = abs(fft(x_recon(:,7,i)));
    plot(f(1:L/2),mag(1:L/2)*2/L,'b--','LineWidth',1.5)
    plot(f(1:L/2),magx(1:L/2)*2/L,'r','LineWidth',1.5)
    xlabel('Frequency(Hz)');ylabel('Mag.');xlim([0 100]);
    legend('Original',['Group',sprintf('%d',i)])
    hold off
end
sgtitle('Floor RF:FFT of SV components')







