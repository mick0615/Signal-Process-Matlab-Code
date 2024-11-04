% 張原嘉-HW1-1-Code
Fs=50;
N1=64;N2=112;N3=128;N=128;
n1=0:N1-1;n2=0:N2-1;n3=0:N3-1;
t1=n1/Fs;t2=n2/Fs;t3=n3/Fs;
ft1=sin(2*pi*3*t1);ft2=sin(2*pi*3*t2);ft3=sin(2*pi*3*t3);
y1=fft(ft1,N);y2=fft(ft2,N);y3=fft(ft3,N);
f=(0:N-1)*Fs/N;
abs1=abs(y1);abs2=abs(y2);abs3=abs(y3);
% N1=64 N2=112 N3=128 plot
figure;
subplot(3,1,1);plot(f(1:N/2),abs1(1:N/2)*2/N);
xlabel('Hz');ylabel('AMP');title('N1=64 N=128');
subplot(3,1,2);plot(f(1:N/2),abs2(1:N/2)*2/N);
xlabel('Hz');ylabel('AMP');title('N1=112 N=128');
subplot(3,1,3);plot(f(1:N/2),abs3(1:N/2)*2/N);
xlabel('Hz');ylabel('AMP');title('N1=N=128');
% NN=512 plot
figure;
NN=512;nn=0:NN-1;tt=nn/Fs;ftt=sin(2*pi*3*tt);
yy=fft(ftt,NN);ff=(0:NN-1)*Fs/NN;abss=abs(yy);
plot(ff(1:NN/2),abss(1:NN/2)*2/NN);
xlabel('Hz');ylabel('AMP');title('N=512');
% N1=64 N2=112 N3=128 NN=512 plot
figure;hold on
plot(f(1:N/2),abs1(1:N/2)*2/N);
xlabel('Hz');ylabel('AMP');title('N1=64 N=128');
plot(f(1:N/2),abs2(1:N/2)*2/N);
xlabel('Hz');ylabel('AMP');title('N1=112 N=128');
plot(f(1:N/2),abs3(1:N/2)*2/N);
xlabel('Hz');ylabel('AMP');title('N1=N=128');
NN=512;nn=0:NN-1;tt=nn/Fs;ftt=sin(2*pi*3*tt);
yy=fft(ftt,NN);ff=(0:NN-1)*Fs/NN;abss=abs(yy);
plot(ff(1:NN/2),abss(1:NN/2)*2/NN);
xlabel('Hz');ylabel('AMP');title('N=64 112 128 512');
legend('N=64','N=112','N=128','N=512');
hold off
% 張原嘉-HW1-2-Code
f1=2;f2=2.1;f3=3.5;
pha1=0;pha2=pi/6;pha3=pi/4;
SR1=50;SR2=10;N1=512;N2=103;
n1=0:N1-1;n2=0:N2-1;
t1=n1/SR1;t2=n2/SR2;
xt1=cos(2*pi*f1*t1+pha1)+1.5*cos(2*pi*f2*t1+pha2)+0.7*cos(2*pi*f3*t1+pha3);
xt2=cos(2*pi*f1*t2+pha1)+1.5*cos(2*pi*f2*t2+pha2)+0.7*cos(2*pi*f3*t2+pha3);
y1=fft(xt1,N1);y2=fft(xt2,N2);
f1=(0:N1-1)*SR1/N1;f2=(0:N2-1)*SR2/N2;
mag1=abs(y1);mag2=abs(y2);
%plot the result
subplot(1,2,1);plot(f1(1:N1/2),mag1(1:N1/2)*2/N1);
title('sampling rate = 50Hz , 512 points');
xlim([0,4.8]);
subplot(1,2,2);plot(f2(1:N2/2),mag2(1:N2/2)*2/N2);
title('sampling rate = 10Hz , 103 points');
xlim([0,4.8]);

