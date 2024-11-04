clc;clear;close all;
% Final Project
%% Load data
data = xlsread('Ming-Li School 1999-9-21 data.xlsx');
input = data(:,28); % accel (gal)
output1 = data(:,8); % accel (gal)
output2 = data(:,15); % accel (gal)
output3 = data(:,21); % accel (gal)
SR = 200;
Fs = 200;
T = 1/ SR;
Time = (0 : 1/SR : (length(data(:,1))-1)*1/SR)';
L = length(Time);
fs = 200; % [Hz], sampling rate
%% high pass and zero mean
%high pass 2Hz Freuency
order = 5;
cff = 2;
[b,a] = butter(order,cff/(SR/2),'low');
inputf = filtfilt(b,a,input);
output1f = filtfilt(b,a,output1);
output2f = filtfilt(b,a,output2);
output3f = filtfilt(b,a,output3);
input = input - inputf;
output1 = output1 - output1f;
output2 = output2 - output2f;
output3 = output3 - output3f;
%zero mean
input = input-mean(input);
output1 = output1-mean(output1);
output2 = output2-mean(output2);
output3 = output3-mean(output3);
%% plot data
% %plot acceleration-time
% figure(1) 
% subplot(3,1,1);plot(Time,output1);grid on
% xlabel('time (sec)');ylabel('Acceleration (gal)');title('ch8 1F樓地板上 川堂');
% subplot(3,1,2);plot(Time,output2);grid on
% xlabel('time (sec)');ylabel('Acceleration (gal)');title('ch15 3F樓地板下 校長室');
% subplot(3,1,3);plot(Time,output3);grid on
% xlabel('time (sec)');ylabel('Acceleration (gal)');title('ch21 RF樓地板下 自然教室');
% saveas(gcf,'Original time data','png')
% sgtitle('acceleration-time')
% 
% %plot FFT
% nFFT = 2^nextpow2(L);
% faxis = (0:nFFT-1)*fs/nFFT;
% FFT1 = abs(fft(output1,nFFT)/(nFFT/2));
% FFT2 = abs(fft(output2,nFFT)/(nFFT/2));
% FFT3 = abs(fft(output3,nFFT)/(nFFT/2));
% figure(2)
% subplot(3,1,1);plot(faxis,FFT1);grid on;xlim([0 25])
% xlabel('frequency(Hz)');ylabel('amplitude');title('ch8 1F樓地板上 川堂')
% subplot(3,1,2);plot(faxis,FFT2);grid on;xlim([0 25])
% xlabel('frequency(Hz)');ylabel('amplitude');title('ch15 3F樓地板下 校長室');
% subplot(3,1,3);plot(faxis,FFT3);grid on;xlim([0 25])
% xlabel('frequency(Hz)');ylabel('amplitude');title('ch21 RF樓地板下 自然教室');
% saveas(gcf,'Original data','png')
% sgtitle('FFT')
%% ARX (0~42sec)
% Take data from 0~42 seconds
indexa = find(Time>=42,1);
% Truncate vectors to that time
inputa = input(1:indexa);
outputa1 = output1(1:indexa);
outputa2 = output2(1:indexa);
outputa3 = output3(1:indexa);

naa = 50; % order of output
nba = 50; % order of input

% Create id data with output and input
outputa = [outputa1 outputa2 outputa3];
Pa = iddata(outputa, inputa, T);

% Perform ARX
sysa1 = arx(Pa(:,1,1),'na', naa, 'nb', nba);
sysa2 = arx(Pa(:,2,1),'na', naa, 'nb', nba);
sysa3 = arx(Pa(:,3,1),'na', naa, 'nb', nba);

% Calculate magnitude and phase
[maga1,phasea1,wa1] = ffplot(sysa1);
[maga2,phasea2,wa2] = ffplot(sysa2);
[maga3,phasea3,wa3] = ffplot(sysa3);

% Create magnitude vector
magnitudea1 = reshape(maga1,1,[]);
magnitudea2 = reshape(maga2,1,[]);
magnitudea3 = reshape(maga3,1,[]);

[a1,b1]=polydata(sysa1);
[a2,b2]=polydata(sysa2);
[a3,b3]=polydata(sysa3);

[r1,p1,k1]=residue(b1,a1);
[r2,p2,k2]=residue(b2,a2);
[r3,p3,k3]=residue(b3,a3);

theta1=angle(p1);
theta2=angle(p2);
theta3=angle(p3);

damp1=log(1./abs(p1))./sqrt(theta1.^2+(log(1./abs(p1))).^2);
damp2=log(1./abs(p2))./sqrt(theta2.^2+(log(1./abs(p2))).^2);
damp3=log(1./abs(p3))./sqrt(theta3.^2+(log(1./abs(p3))).^2);

freq1=log(1./abs(p1))./(2*pi*damp1*T);
freq2=log(1./abs(p2))./(2*pi*damp2*T);
freq3=log(1./abs(p3))./(2*pi*damp3*T);

[freq1,index1]=sort(freq1);
[freq2,index2]=sort(freq2);
[freq3,index3]=sort(freq3);

for i=1:length(index1)
    damping1(i)=damp1(index1(i));
    damping2(i)=damp2(index2(i));
    damping3(i)=damp3(index3(i));
end
[pksa1,pksa1id] = findpeaks(magnitudea1);
[pksa2,pksa2id] = findpeaks(magnitudea2);
[pksa3,pksa3id] = findpeaks(magnitudea3);

% Plot FRF
figure;
semilogy(wa1,magnitudea1)
hold on
plot(wa1(pksa1id(1:3)),magnitudea1(pksa1id(1:3)),"ro")
text(wa1(pksa1id(1)),magnitudea1(pksa1id(1)),[sprintf(' %.2f', wa1(pksa1id(1))),'Hz'])
text(wa1(pksa1id(2)),magnitudea1(pksa1id(2)),[sprintf(' %.2f', wa1(pksa1id(2))),'Hz'])
text(wa1(pksa1id(3)),magnitudea1(pksa1id(3)),[sprintf(' %.2f', wa1(pksa1id(3))),'Hz'])
hold off
grid on; grid minor
xlabel('Frequency [Hz]')
ylabel('Amplitude (A.U.)')
%xlim([0 19])
% ylim([1e-8 1e+3])
title('Frequency Response Function ch8 1F樓地板上 川堂 0~42 sec')

figure;
semilogy(wa2,magnitudea2)
hold on
plot(wa2(pksa2id(1:3)),magnitudea2(pksa2id(1:3)),"ro")
text(wa2(pksa2id(1)),magnitudea2(pksa2id(1)),[sprintf(' %.2f', wa1(pksa2id(1))),'Hz'])
text(wa2(pksa2id(2)),magnitudea2(pksa2id(2)),[sprintf(' %.2f', wa1(pksa2id(2))),'Hz'])
text(wa2(pksa2id(3)),magnitudea2(pksa2id(3)),[sprintf(' %.2f', wa1(pksa2id(3))),'Hz'])
hold off
grid on; grid minor
xlabel('Frequency [Hz]')
ylabel('Amplitude (A.U.)')
%xlim([0 19])
% ylim([1e-8 1e+3])
title('Frequency Response Function ch15 3F樓地板下 校長室 0~42 sec')
figure;
semilogy(wa3,magnitudea3)
hold on
plot(wa3(pksa3id(1:3)),magnitudea3(pksa3id(1:3)),"ro")
text(wa3(pksa3id(1)),magnitudea3(pksa3id(1)),[sprintf(' %.2f', wa3(pksa3id(1))),'Hz'])
text(wa3(pksa3id(2)),magnitudea3(pksa3id(2)),[sprintf(' %.2f', wa3(pksa3id(2))),'Hz'])
text(wa3(pksa3id(3)),magnitudea3(pksa3id(3)),[sprintf(' %.2f', wa3(pksa3id(3))),'Hz'])
hold off
grid on; grid minor
xlabel('Frequency [Hz]')
ylabel('Amplitude (A.U.)')
%xlim([0 19])
% ylim([1e-8 1e+3])
title('Frequency Response Function ch21 RF樓地板下 自然教室 0~42 sec')

% Plot transfer functions
figure;
semilogy(wa1 , magnitudea1,'LineWidth',1) ; hold on
[~,locs] = findpeaks(magnitudea1);
legendInfo{1} = ['ch8 1F樓地板上 川堂  /  ',num2str(round(wa1(locs(1)),2)),'Hz  /  ',num2str(round(damping1(1)*100,2))];
[~,seat1] = min(abs(wa1-4)); 
temp_Axis(1,1) = min(magnitudea1(1:seat1+1));
temp_Axis(2,1) = max(magnitudea1(1:seat1+1));
xlabel('Freq. (Hz)') ; ylabel('Amplitude') ; grid on
axis([0 4 min(temp_Axis(1,:))*0.5 max(temp_Axis(2,:))*5])
title(['Frequency Response Function from ARX (50,50)'])
plot([0,0],'color','none') ; 
legendInfo{end+1} = '樓層 / 自然頻率 / 阻尼比(%)';
legend(legendInfo)
set(gcf,'unit','normalized','position',[0.3,0.3,0.3,0.4]); set(gca,'fontsize',11)
hold off
figure;
semilogy(wa2 , magnitudea2,'LineWidth',1) ; hold on
[~,locs] = findpeaks(magnitudea2);
legendInfo{1} = ['ch15 3F樓地板下 校長室  /  ',num2str(round(wa2(locs(1)),2)),'Hz  /  ',num2str(round(damping2(1)*100,2))];
[~,seat1] = min(abs(wa2-4)); 
temp_Axis(1,1) = min(magnitudea2(1:seat1+1));
temp_Axis(2,1) = max(magnitudea2(1:seat1+1));
xlabel('Freq. (Hz)') ; ylabel('Amplitude') ; grid on
axis([0 4 min(temp_Axis(1,:))*0.5 max(temp_Axis(2,:))*5])
title(['Frequency Response Function from ARX (50,50)'])
plot([0,0],'color','none') ; 
legendInfo{end+1} = '樓層 / 自然頻率 / 阻尼比(%)';
legend(legendInfo)
set(gcf,'unit','normalized','position',[0.3,0.3,0.3,0.4]); set(gca,'fontsize',11)
hold off
figure;
semilogy(wa3 , magnitudea3,'LineWidth',1) ; hold on
[~,locs] = findpeaks(magnitudea3);
legendInfo{1} = ['ch21 RF樓地板下 自然教室  /  ',num2str(round(wa3(locs(1)),2)),'Hz  /  ',num2str(round(damping3(1)*100,2))];
[~,seat1] = min(abs(wa3-4)); 
temp_Axis(1,1) = min(magnitudea3(1:seat1+1));
temp_Axis(2,1) = max(magnitudea3(1:seat1+1));
xlabel('Freq. (Hz)') ; ylabel('Amplitude') ; grid on
axis([0 4 min(temp_Axis(1,:))*0.5 max(temp_Axis(2,:))*5])
title(['Frequency Response Function from ARX (50,50)'])
plot([0,0],'color','none') ; 
legendInfo{end+1} = '樓層 / 自然頻率 / 阻尼比(%)';
legend(legendInfo)
set(gcf,'unit','normalized','position',[0.3,0.3,0.3,0.4]); set(gca,'fontsize',11)
hold off

%% ARX (0~102sec)
% Take data from 0~102 seconds

% Truncate vectors to that time
inputb = input(1:end);
outputb1 = output1(1:end);
outputb2 = output2(1:end);
outputb3 = output3(1:end);

nab = 50; % order of output (user defined parameter)
nbb = 50; % order of input (user defined parameter)

% Create id data with output and input
outputb = [outputb1,outputb2,outputb3];
Pb = iddata(outputb, inputb, T);

% Perform ARX
sysb1 = arx(Pb(:,1,1),'na', nab, 'nb', nbb);
sysb2 = arx(Pb(:,2,1),'na', nab, 'nb', nbb);
sysb3 = arx(Pb(:,3,1),'na', nab, 'nb', nbb);

% % Return polynomial data
% [ab,bb] = ploydata(sysb);

% Calculate magnitude and phase
[magb1,phaseb1,wb1] = ffplot(sysb1);
[magb2,phaseb2,wb2] = ffplot(sysb2);
[magb3,phaseb3,wb3] = ffplot(sysb3);

% Create magnitude vector
magnitudeb1 = reshape(magb1,1,[]);
magnitudeb2 = reshape(magb2,1,[]);
magnitudeb3 = reshape(magb3,1,[]);

[a1,b1]=polydata(sysb1);
[a2,b2]=polydata(sysb2);
[a3,b3]=polydata(sysb3);

[r1,p1,k1]=residue(b1,a1);
[r2,p2,k2]=residue(b2,a2);
[r3,p3,k3]=residue(b3,a3);

theta1=angle(p1);
theta2=angle(p2);
theta3=angle(p3);

damp1=log(1./abs(p1))./sqrt(theta1.^2+(log(1./abs(p1))).^2);
damp2=log(1./abs(p2))./sqrt(theta2.^2+(log(1./abs(p2))).^2);
damp3=log(1./abs(p3))./sqrt(theta3.^2+(log(1./abs(p3))).^2);

freq1=log(1./abs(p1))./(2*pi*damp1*T);
freq2=log(1./abs(p2))./(2*pi*damp2*T);
freq3=log(1./abs(p3))./(2*pi*damp3*T);

[freq1,index1]=sort(freq1);
[freq2,index2]=sort(freq2);
[freq3,index3]=sort(freq3);

for i=1:length(index1)
    damping1(i)=damp1(index1(i));
    damping2(i)=damp2(index2(i));
    damping3(i)=damp3(index3(i));
end
[pksb1,pksb1id] = findpeaks(magnitudeb1);
[pksb2,pksb2id] = findpeaks(magnitudeb2);
[pksb3,pksb3id] = findpeaks(magnitudeb3);
% % Partial Fraction Decomposition
% [rb,pb,kb] = residue(bb,ab);

% Plot FRF
% Plot results
figure;
semilogy(wb1,magnitudeb1)
hold on
plot(wb1(pksb1id(1:3)),magnitudeb1(pksb1id(1:3)),"ro")
text(wb1(pksb1id(1)),magnitudeb1(pksb1id(1)),[sprintf(' %.2f', wb1(pksb1id(1))),'Hz'])
text(wb1(pksb1id(2)),magnitudeb1(pksb1id(2)),[sprintf(' %.2f', wb1(pksb1id(2))),'Hz'])
text(wb1(pksb1id(3)),magnitudeb1(pksb1id(3)),[sprintf(' %.2f', wb1(pksb1id(3))),'Hz'])
hold off
grid on; grid minor
xlabel('Frequency [Hz]')
ylabel('Amplitude (A.U.)')
%xlim([0 19])
% ylim([1e-8 1e+3])
title('Frequency Response Function ch8 1F樓地板上 川堂 0~102 sec')

% Plot results
figure;
semilogy(wb2,magnitudeb2)
hold on
plot(wb2(pksb2id(1:3)),magnitudeb2(pksb2id(1:3)),"ro")
text(wb2(pksb2id(1)),magnitudeb2(pksb2id(1)),[sprintf(' %.2f', wb2(pksb2id(1))),'Hz'])
text(wb2(pksb2id(2)),magnitudeb2(pksb2id(2)),[sprintf(' %.2f', wb2(pksb2id(2))),'Hz'])
text(wb2(pksb2id(3)),magnitudeb2(pksb2id(3)),[sprintf(' %.2f', wb2(pksb2id(3))),'Hz'])
hold off
grid on; grid minor
xlabel('Frequency [Hz]')
ylabel('Amplitude (A.U.)')
%xlim([0 19])
% ylim([1e-8 1e+3])
title('Frequency Response Function ch15 3F樓地板下 校長室 0~102 sec')
% Plot results
figure;
semilogy(wb3,magnitudeb3)
hold on
plot(wb3(pksb3id(1:3)),magnitudeb3(pksb3id(1:3)),"ro")
text(wb3(pksb3id(1)),magnitudeb3(pksb3id(1)),[sprintf(' %.2f', wb3(pksb3id(1))),'Hz'])
text(wb3(pksb3id(2)),magnitudeb3(pksb3id(2)),[sprintf(' %.2f', wb3(pksb3id(2))),'Hz'])
text(wb3(pksb3id(3)),magnitudeb3(pksb3id(3)),[sprintf(' %.2f', wb3(pksb3id(3))),'Hz'])
hold off
grid on; grid minor
xlabel('Frequency [Hz]')
ylabel('Amplitude (A.U.)')
%xlim([0 19])
% ylim([1e-8 1e+3])
title('Frequency Response Function ch21 RF樓地板下 自然教室 0~102 sec')

% Plot transfer functions
figure;
semilogy(wb1 , magnitudeb1,'LineWidth',1) ;hold on
[~,locs] = findpeaks(magnitudeb1);
legendInfo{1} = ['ch8 1F樓地板上 川堂  /  ',num2str(round(wb1(locs(1)),2)),'Hz  /  ',num2str(round(damping1(1)*100,2))];
[~,seat1] = min(abs(wb1-4)); 
temp_Axis(1,1) = min(magnitudeb1(1:seat1+1));
temp_Axis(2,1) = max(magnitudeb1(1:seat1+1));
xlabel('Freq. (Hz)') ; ylabel('Amplitude') ; grid on
axis([0 4 min(temp_Axis(1,:))*0.5 max(temp_Axis(2,:))*5])
title(['Frequency Response Function from ARX (50,50)'])
plot([0,0],'color','none') ; 
legendInfo{end+1} = '樓層 / 自然頻率 / 阻尼比(%)';
legend(legendInfo)
set(gcf,'unit','normalized','position',[0.3,0.3,0.3,0.4]); set(gca,'fontsize',11)
figure;
semilogy(wb2 , magnitudeb2,'LineWidth',1) ;hold on
[~,locs] = findpeaks(magnitudeb2);
legendInfo{1} = ['ch15 3F樓地板下 校長室  /  ',num2str(round(wb2(locs(1)),2)),'Hz  /  ',num2str(round(damping2(1)*100,2))];
[~,seat1] = min(abs(wb2-4)); 
temp_Axis(1,1) = min(magnitudeb2(1:seat1+1));
temp_Axis(2,1) = max(magnitudeb2(1:seat1+1));
xlabel('Freq. (Hz)') ; ylabel('Amplitude') ; grid on
axis([0 4 min(temp_Axis(1,:))*0.5 max(temp_Axis(2,:))*5])
title(['Frequency Response Function from ARX (50,50)'])
plot([0,0],'color','none') ; 
legendInfo{end+1} = '樓層 / 自然頻率 / 阻尼比(%)';
legend(legendInfo)
set(gcf,'unit','normalized','position',[0.3,0.3,0.3,0.4]); set(gca,'fontsize',11)
figure;
semilogy(wb3 , magnitudeb3,'LineWidth',1) ;hold on
[~,locs] = findpeaks(magnitudeb3);
legendInfo{1} = ['ch21 RF樓地板下 自然教室  /  ',num2str(round(wb3(locs(1)),2)),'Hz  /  ',num2str(round(damping3(1)*100,2))];
[~,seat1] = min(abs(wb3-4)); 
temp_Axis(1,1) = min(magnitudeb3(1:seat1+1));
temp_Axis(2,1) = max(magnitudeb3(1:seat1+1));
xlabel('Freq. (Hz)') ; ylabel('Amplitude') ; grid on
axis([0 4 min(temp_Axis(1,:))*0.5 max(temp_Axis(2,:))*5])
title(['Frequency Response Function from ARX (50,50)'])
plot([0,0],'color','none') ; 
legendInfo{end+1} = '樓層 / 自然頻率 / 阻尼比(%)';
legend(legendInfo)
set(gcf,'unit','normalized','position',[0.3,0.3,0.3,0.4]); set(gca,'fontsize',11)
%% mode shape(0~102sec)
for i = 1:3
mode_story1(i) = magnitudeb1(pksb1id(i)).*exp(1i*phaseb1(pksb1id(i))) ;
mode_story2(i) = magnitudeb2(pksb2id(i)).*exp(1i*phaseb2(pksb2id(i))) ;
mode_story3(i) = magnitudeb3(pksb3id(i)).*exp(1i*phaseb3(pksb3id(i))) ;
end

mode1 = [0,mode_story1(1),mode_story2(1),mode_story3(1)] ;
mode2 = [0,mode_story1(2),mode_story2(2),mode_story3(2)] ;
mode3 = [0,mode_story1(3),mode_story2(3),mode_story3(3)] ;
vertical = zeros(4,1);
floorplot = (0:1:3);

%第一模態不考慮延遲問題 直接從Frequency Response Function 求 Mode Shape
figure;
fr1 = [0;3.55;5.489;9.16];
fr1 = [0;fr1(2)/fr1(4);fr1(3)/fr1(4);fr1(4)/fr1(4)];
plot(vertical, floorplot, '--r', fr1, floorplot, '-bo')
title('Mode Shape 1 (0~102sec)');xlim([-3 3]);

figure;
for i = 1:4
    fr2(i) = sign(real(mode2(:,i)))*abs(mode2(:,i)); 
end
fr2 = [0;fr2(2)/fr2(4);fr2(3)/fr2(4);fr2(4)/fr2(4)];
plot(vertical, floorplot, '--r', fr2, floorplot, '-bo')
title('Mode Shape 2 (0~102sec)');xlim([-15 15]);

figure;
for i = 1:4
    fr3(i) = sign(real(mode3(:,i)))*abs(mode3(:,i)); 
end
fr3 = [0;fr3(2)/fr3(4);fr3(3)/fr3(4);fr3(4)/fr3(4)];
plot(vertical, floorplot, '--r', fr3, floorplot, '-bo')
title('Mode Shape 3 (0~102sec)');xlim([-8 8]);