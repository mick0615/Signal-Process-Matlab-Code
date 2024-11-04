% Load Data
clc;clear;close all;
data = xlsread('Ming-Li School 1999-9-21 data.xlsx');
SR = 200 ;
time = (0 : 1/SR : (length(data(:,1))-1)*1/SR)';
input = data(:,28); % accel (gal)
output = data(:,21); % accel (gal)
fs = 200; % [Hz], sampling rate
T = 1/fs; % [seconds], sampling period
L = length(time); % length of signal

% %low pass 20Hz Freuency
% order = 6;
% cff = 20;
% [b,a] = butter(order,cff/(SR/2),'low');
% input = filtfilt(b,a,input);
% output = filtfilt(b,a,output);
% %zero mean
% input = input-mean(input);
% output = output-mean(output);


%% MCMW-VCF
% FreqAnal = (0:0.01:20)';
% sigma = 5.0;
% type_coef = 'energy';
% [coefs1,psi,T1]=MCMW_VCF(output,FreqAnal,fs,type_coef,2,time);
% figure(3) %MCMW-VCF Story 1
% imagesc(T1,FreqAnal,coefs1);grid on;set(gca,'YDir','normal');colormap jet;
% xlabel('Time(sec)');ylabel('Freq. (Hz)')
% title("MCMW-VCF ( Story 1, \sigma = "+sigma+" )")
% ylim([0 10])
% set(gcf,'unit','normalized','position',[0.2,0.2,0.2,0.3]);
% colorbar;
% %% ARX 0 ~ 42 sec
% % Take data from 0~42 seconds
% indexa = find(time>=42,1);
% 
% % Truncate vectors to that time
% inputa = input(1:indexa);
% outputa = output(1:indexa);
% naa = 30; % order of output
% nba = 30; % order of input
% 
% % Create id data with output and input
% Pa = iddata(outputa, inputa, T);
% 
% % Perform ARX
% sysa = arx(Pa,'na', naa, 'nb', nba);
% 
% % % Return polynomial data
% %[aa,ba] = ploydata(sysa);
% 
% % Calculate magnitude and phase
% [maga1,phasea1,wa1] = ffplot(sysa);
% 
% % Create magnitude vector
% magnitudea1 = reshape(maga1,1,[]);
% [pksa1,pksa1id] = findpeaks(magnitudea1);
% 
% % % Partial Fraction Decomposition
% %[ra,pa,ka] = residue(ba,aa);
% 
% % Plot results
% figure;
% semilogy(wa1,magnitudea1)
% hold on
% plot(wa1(pksa1id(1:3)),magnitudea1(pksa1id(1:3)),"ro")
% text(wa1(pksa1id(1)),magnitudea1(pksa1id(1)),[sprintf(' %.2f', wa1(pksa1id(1))),'Hz'])
% text(wa1(pksa1id(2)),magnitudea1(pksa1id(2)),[sprintf(' %.2f', wa1(pksa1id(2))),'Hz'])
% text(wa1(pksa1id(3)),magnitudea1(pksa1id(3)),[sprintf(' %.2f', wa1(pksa1id(3))),'Hz'])
% hold off
% grid on; grid minor
% xlabel('Frequency [Hz]')
% ylabel('Amplitude (A.U.)')
% %xlim([0 19])
% % ylim([1e-8 1e+3])
% title('Frequency Response Function 1999 Chichi Earthquake (0~42 seconds)')
% legend('ch21 RF樓地板下 自然教室','location','northeast')
% %% ARX 42sec ~ end
% % Take data from 42~end seconds
% % Truncate vectors to that time
% inputb = input(indexa:end);
% outputb = output(indexa:end);
% 
% nab = 30; % order of output
% nbb = 30; % order of input
% 
% % Create id data with output and input
% Pb = iddata(outputb, inputb, T);
% 
% % Perform ARX
% sysb1 = arx(Pb,'na', nab, 'nb', nbb);
% 
% % % Return polynomial data
% %[ab,bb] = ploydata(sysb);
% 
% % Calculate magnitude and phase
% [magb1,phaseb1,wb1] = ffplot(sysb1);
% 
% % Create magnitude vector
% magnitudeb1 = reshape(magb1,1,[]);
% [pksb1,pksb1id] = findpeaks(magnitudeb1);
% 
% % % Partial Fraction Decomposition
% % [rb,pb,kb] = residue(bb,ab);
% 
% % Plot results
% figure;
% semilogy(wb1,magnitudeb1)
% hold on
% plot(wb1(pksb1id(1:3)),magnitudeb1(pksb1id(1:3)),"ro")
% text(wb1(pksb1id(1)),magnitudeb1(pksb1id(1)),[sprintf(' %.2f', wb1(pksb1id(1))),'Hz'])
% text(wb1(pksb1id(2)),magnitudeb1(pksb1id(2)),[sprintf(' %.2f', wb1(pksb1id(2))),'Hz'])
% text(wb1(pksb1id(3)),magnitudeb1(pksb1id(3)),[sprintf(' %.2f', wb1(pksb1id(3))),'Hz'])
% hold off
% grid on; grid minor
% xlabel('Frequency [Hz]')
% ylabel('Amplitude (A.U.)')
% %xlim([0 19])
% % ylim([1e-8 1e+3])
% title('Frequency Response Function 1999 Chichi Earthquake (42~102 seconds)')
% legend('ch21 RF樓地板下 自然教室','location','northeast')

%% PROBLEM 3: Recursive ARS model (RARX)
% RARX approach 1 ===================
% RecursiveARX model
% indexa = find(time>=42,1);
% indexc = indexa;
% indexd = find(time>=30);
% output3d2 = output(indexc:indexd);
% inputd2 = input(indexc:indexd);
indexc = 1;
indexd = 20480;
nn = [16 16 0];
obj = recursiveARX(nn);
A = zeros(length(time),nn(1)+1);
B = zeros(length(time),nn(2));
output = data(:,21);
 for itime = 1:1:length(time)
      [A(itime,:),B(itime,:),out] = step(obj,output(itime),input(itime));
 end


figure;
plot(time(indexc:indexd),A(indexc:indexd,:))
grid on; grid minor
xlim([9.5 20])
title('Coefficients (a0~a17) calculated by ARX model')
ylabel('Coefficient Value')
xlabel('Time [sec]')

figure;
plot(time(indexc:indexd),B(indexc:indexd,:))
grid on; grid minor
xlim([9.5 20])
title('Coefficients (b0~b17) calculated by ARX model')
ylabel('Coefficient Value')
xlabel('Time [sec]')

% Reconstruct Transfor Function
% Create Plolynomial at selected time values
% Construct the input/output model
sysd31 = idpoly(A(fs*30.0+1,:),B(fs*30.0+1,:),[],[],[],[],T); % 50 Hz * 8.50 sec = 300 pts
sysd32 = idpoly(A(fs*40.0+1,:),B(fs*40.0+1,:),[],[],[],[],T); % 50 Hz * 9.56 sec = 400 pts
sysd33 = idpoly(A(fs*60.0+1,:),B(fs*60.0+1,:),[],[],[],[],T); % 50 Hz * 20.0 sec = 300 pts
sysd34 = idpoly(A(fs*70.0+1,:),B(fs*70.0+1,:),[],[],[],[],T); % 50 Hz * 30.0 sec = 400 pts

% Calculated magnitude and phase for each
[magd1,phased1,wd1] = ffplot(sysd31);
[magd2,phased2,wd2] = ffplot(sysd32);
[magd3,phased3,wd3] = ffplot(sysd33);
[magd4,phased4,wd4] = ffplot(sysd34);

% Create magnitude vector for each
magnituded1 = reshape(magd1,1,[]);
magnituded2 = reshape(magd2,1,[]);
magnituded3 = reshape(magd3,1,[]);
magnituded4 = reshape(magd4,1,[]);

[pksd1,pksd1id] = findpeaks(magnituded1);
[pksd2,pksd2id] = findpeaks(magnituded2);
[pksd3,pksd3id] = findpeaks(magnituded3);
[pksd4,pksd4id] = findpeaks(magnituded4);

pkwd1 = wd1(pksd1id);
pkwd2 = wd2(pksd2id);
pkwd3 = wd3(pksd3id);
pkwd4 = wd4(pksd4id);

% Plot results
figure;
semilogy(wd1,magnituded1)
hold on
semilogy(wd2,magnituded2,"--")
semilogy(wd3,magnituded3)
semilogy(wd4,magnituded4,"--")
hold off
grid on; grid minor
xlabel('Frequency ([Hz]')
ylabel('Amplitude')
xlim([0 30])
title('Frequency Response 1999 Chichi Earthquake at Selected Instant of Times')
legend("at t = 30.0 sec","at t = 40.0 sec","at t = 60.0 sec","at t = 70.0 sec")
%% RARX approach 2 ===================
fs = 200;
x = input; % input inputb
y = output; % output story3b
z = [y x]; % [output, input]
forgetFactor = 1;
% Examine parameter a and b
na = 30;
nb = 30;
[thm, noise] = rarx(z, [na,nb,1],'ff',forgetFactor);
a = thm(:,1:na);
b = thm(:,na+1:na+nb);
% timed = time(indexc:indexd);
timed = time;

figure()
plot(timed,a)
title('Parameters a')
xlabel('T [sec]')
ylabel('a');ylim([-2.5 2.5]);xlim([0 102.39])
grid on; grid minor

figure()
plot(timed,b)
title('Parameters b')
xlabel('T [sec]')
ylabel('b');ylim([-0.5 0.5]);xlim([0 102.39])
grid on; grid minor
SR = fs;
%%% Find the modal frequencies -------------------
na = 30;
nb = 30;
num = [na,nb,1];
[thm,noise] = rarx(z,num,'ff',forgetFactor);
n = 100;
for itime = (n+num(1)+num(2)):length(timed)
    a = [1,thm(itime,1:num(1))];
    b = [0,thm(itime,(num(1)+1):end)];
    p = roots(a);
    kesi = log(p)*fs;
    norm = -real(kesi)./abs(kesi);
    [tmp,tmp_index] = sort(abs(kesi)/(2*pi),'ascend');
    fk(itime,:) = tmp; % natural freq
    tmp = norm(tmp_index); % damping ratio
    save_a(:,:,itime) = a;
    save_b(:,:,itime) = b;
    save_p(:,:,itime) = p;
end
freq = [fk(:,1) fk(:,3) fk(:,5)];
index = find(freq(:,1)==0);
ms_index = freq(size(index,1)+1,:);
ms1 = zeros(size(index,1),3);
ms2 = freq(size(index,1)+1:end,:);
for itime = 1:size(ms1,1)
    ms1(itime,:) = ms_index;
end
ms = [ms1;ms2];
% plot - Time-Varying Modal Frequencie  -------------------
figure;
grid on; grid minor;
hold on
set(gcf,'color','w')
plot(timed,ms)
legend('Mode 1','Mode 2','Mode 3')
xlabel('Time [s]');
ylabel('Modal Frequencies [Hz]')
xlim([0 time(end)]);
%ylim([0 6])
title("Time-Varying Modal Frequencies of the System (Forgeting Factor = "+forgetFactor+" )")
%%% Plot - System poles  -------------------
figure()
[hz1, hp1, ht1] = zplane(save_p(:,:,4000)); % 200Hz*20sec
hold on
[hz2, hp2, ht2] = zplane(save_p(:,:,12000)); % 200Hz*60sec
hold off
set(findobj(hz1, 'Type', 'line'), 'Color', 'b','displayname',"time = 20 sec");
set(findobj(hp1, 'Type', 'line'),'HandleVisibility','off');
set(findobj(ht1, 'Type', 'line'),'HandleVisibility','off');
set(findobj(hz2, 'Type', 'line'), 'Color', 'r','displayname',"time = 60 sec");
set(findobj(hp2, 'Type', 'line'),'HandleVisibility','off');
set(findobj(ht2, 'Type', 'line'),'HandleVisibility','off');
xlim([-1.3 1.3])
ylim([-1.3 1.3])
title("System Poles")
grid on; grid minor
legend