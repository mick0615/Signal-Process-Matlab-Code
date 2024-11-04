% M11205314-HW1-2-Code
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

