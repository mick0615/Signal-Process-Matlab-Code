clc , clear, close all

% Use relative acceleration to estimate displacement with permanent
% deformation
data = xlsread('Ming-Li School 1999-9-21 data.xlsx');
SR = 200 ;
Time = (0 : 1/SR : (length(data(:,1))-1)*1/SR)';
base_acc = data(:,28);
roof_acc = data(:,21);
%Reladis=p1(:,5);

Accel = roof_acc - base_acc ;

%FigNo = 1 ;
%mkdir(fullfile(cd,['ResultOutput','\',StationName])) ; SavePath = fullfile(cd,['ResultOutput','\',StationName]);
SavePath = cd; % 新增******

figure (1)
plot( Time , Accel , 'LineWidth', 1)
xlabel('Time (sec)') ; ylabel(' Accel') ; 
xlim([Time(1) Time(end)])
set(gca,'xminortick','on') ; set(gca , 'FontSize',12) ; title('Original-Data','FontSize',14)
set(gcf,'unit','normalized','position',[0.3,0.3,0.4,0.3]);
%saveas(gcf,fullfile(cd,'Ori-Data.png'))


% -------------------
% 永久變位分析 (building structural response)
% -------------------
% [1] SSA分析 - 獲得重組訊號

L = 2048;
AmpAccel = Accel + 10*max(abs(Accel));
[RecAccel] = SSA(SR , AmpAccel , size(Accel,1) , size(Accel,2) , L , size(Accel,1)-L+1 , SavePath);

%check 重組訊號與原訊號是否相似
figure (2)
plot(Time , Accel); hold on ; plot(Time , RecAccel);
legend('Comparison-Data' , '重組資料','fontsize',12)
title(' check 重組訊號與原訊號是否相似','fontsize',14)

% RecAccel=Accel;

% [2-1] EMD
[imf,~,~] = emd(RecAccel,'Interpolation','pchip');

figure (3)
for j = 1:size(imf,2)
    subplot(size(imf,2),1,j)
    plot(Time,imf(:,j)); grid on;
    ylabel(['IMF',num2str(j)]) ; xlim([Time(1) Time(end)])
    if j ==1
        title('EMD(acc)','FontSize',14)
    end
end
xlabel('Time (sec)'); set(gcf,'unit','normalized','position',[0.1,0.1,0.4,0.75]);
% saveas(gcf,fullfile(SavePath,['Fig5_EMD.png'])) ; FigNo = FigNo + 1;

figure (4)
for j = 1:size(imf,2)
    n = length(imf(:,j));
    f = 0:(SR/n):(SR/2-SR/n);
    y = abs(fft(imf(:,j))/n);
    y = y(1:n/2);
    subplot(size(imf,2),1,j);
    plot(f,y); grid on;
    xlim([0 20]); ylabel(['IMF',num2str(j)])
    if j ==1
        title('EMD FFT(acc)','FontSize',14) ;
    end
end
xlabel('Frequency (Hz)'); set(gcf,'unit','normalized','position',[0.1,0.1,0.4,0.75]); clear n f y
% saveas(gcf,fullfile(SavePath,['Fig6_EMD_FFT.png'])) 
clear L n f y faxis nfft N AxisValue j Signal_FFT

% [2-2] Calculate Energy To check A1
% 計算energy是為了判斷Acc要取到第幾個IMF(從圖來判斷IMF到第幾個，能量沒有明顯增加)
for j = 1 : size(imf,2)
    a(:,j) = sum(imf(:,1:j),2);
    E(:,j) = cumtrapz(Time,RecAccel.*a(:,j));
end
figure (5)
for j = 1 : size(E,2)
    if j <= 7
        plot(Time , E(:,j) , 'LineWidth',1) ; grid on ; hold on ;
    else
        plot(Time , E(:,j) , '--','LineWidth',1);
    end
end
xlabel('Time (sec)'); ylabel('Energy'); xlim([Time(1) Time(end)])
legend('IMF1','IMF1+IMF2','IMF1+...+IMF3','IMF1+...+IMF4','IMF1+...+IMF5','IMF1+...+IMF6','IMF1+...+IMF7','IMF1+...+IMF8','IMF1+...+IMF9','IMF1+...+IMF10','Location','northeastoutside','fontsize',11)
set(gcf,'unit','normalized','position',[0.3,0.3,0.35,0.45]);
% saveas(gcf,fullfile(SavePath,['Fig7_Energy.png']))
clear a E

% [2-2-L] A1
N = input("累積能量A1累積至IMF_N: N = "); % 由"Energy"決定
A1 = sum(imf(:,1:N),2);
figure (6)
subplot(2,1,1)
plot(Time , A1) ; grid on ;
ylabel('Acceleration (gal)') ; xlabel('Time (sec)') ;
title(['A_1 = IMF1+...+IMF',num2str(N)],'FontSize',14) ;
xlim([Time(1) Time(end)]) ;

n = length(A1);
f = 0:(SR/n):(SR/2-SR/n);
y = abs(fft(A1)/n);
y = y(1:n/2);
subplot(2,1,2)
plot(f,y); grid on;
xlabel('Frequency (Hz)'); ylabel('Amplitude') ;
title('A_1 FFT ', 'FontSize',14) ; xlim([0 20]);
set(gcf,'unit','normalized','position',[0.3,0.3,0.35,0.45]);
% saveas(gcf,fullfile(SavePath,['Fig8_A1_Accel&FFT.png'])) 
clear n y f

fcut = input("fcut = "); % 由"A1傅氏譜"決定
ft = fcut + 0.1;
[D1] = drift_ormsby(A1 , SR , Time , fcut , ft); % 位移，單位mm

DT=D1;
D1=detrend(DT, 6);

figure (7)
plot(Time , D1 , 'LineWidth',1); grid on ; hold on ;
% plot(Time , LVDT(:,i) ,'r-' , 'LineWidth' ,1);  %--------------------- 沒有資料*
ylabel('Displacement (mm)'); xlabel('Time (sec)');
axis([Time(1) Time(end) max(abs(D1))*-1.1  max(abs(D1))*1.1]) ; legend('D_1','LVDT');
title('Displacement  D_1','fontsize',14)
set(gcf,'unit','normalized','position',[0.2,0.2,0.4,0.25]);
% saveas(gcf,fullfile(SavePath,['Fig9_D1.png'])) 
clear N ;

% [2-2-R] A2
[y] = ormsby(A1 , ft , fcut , SR , length(A1) );
A1R = A1-y; clear y
A2 = A1R + (RecAccel - A1); % A2 = A1R + A2;
figure (8)
plot(Time , A2); grid on;
xlabel('Time (sec)'); ylabel('Acceleration (gal)') ;
xlim([Time(1) Time(end)]) ; title('A_2 ','FontSize',14);
set(gcf,'unit','normalized','position',[0.2,0.2,0.4,0.25]);
% saveas(gcf,fullfile(SavePath,['Fig10_A2.png']))

figure (9)
subplot(3,1,1)
plot(Time , A2); grid on;
ylabel('Acceleration (gal)'); set(gca,'fontsize',10)
axis([Time(1) Time(end) max(abs(A2))*-1.1  max(abs(A2))*1.1]) ;
subplot(3,1,2)
plot(Time , cumtrapz(Time , A2)); grid on;
ylabel('Velocity (cm/sec)'); set(gca,'fontsize',10)
axis([Time(1) Time(end) max(abs(cumtrapz(Time , A2)))*-1.1  max(abs(cumtrapz(Time , A2)))*1.1]) ;
subplot(3,1,3);
plot(Time , cumtrapz(Time , cumtrapz(Time,A2))*10) ; grid on ; hold on;
axis([Time(1) Time(end) max(abs(cumtrapz(Time , cumtrapz(Time,A2))*10))*-1.1  max(abs(cumtrapz(Time , cumtrapz(Time,A2))*10))*1.1]) ;
set(gca,'fontsize',10)
xlabel('Time (sec)'); ylabel('Displacement (mm)'); sgtitle('\bfA_2','FontSize',14);
set(gcf,'unit','normalized','position',[0.3,0.3,0.35,0.45]);
% saveas(gcf,fullfile(SavePath,['Fig11_A2_AVD.png']))

% [4-1-R] v2 - EMD
V2 = cumtrapz(Time , A2);
[imf2 , r , ~] = emd(V2,'Interpolation','pchip');

figure (10)
for j = 1 : size(imf2,2)
    subplot(size(imf2,2),1,j)
    plot(Time , imf2(:,j) , 'LineWidth',1); grid on;
    ylabel(['IMF',num2str(j)]); xlim([Time(1) Time(end)])
end
xlabel('Time (sec)'); set(gcf,'unit','normalized','position',[0.1,0.1,0.4,0.75]); sgtitle('\bfṼ_2','fontsize',14);
% saveas(gcf,fullfile(SavePath,['Fig12_V2_EMD.png'])) 

for j = 1 : size(imf2,2)
    vel_accumu(:,j) = sum(imf2(:,1:j),2);
end
yaxisname = {'IMF1','IMF1+IMF2','IMF1+...+IMF3','IMF1+...+IMF4','IMF1+...+IMF5','IMF1+...+IMF6','IMF1+...+IMF7','IMF1+...+IMF8'};

figure (11)
for j = 1 : size(vel_accumu,2)
    subplot(size(vel_accumu,2),1,j)
    plot(Time,vel_accumu(:,j),LineWidth=1); grid on
    ylabel(yaxisname{j});xlim([Time(1) Time(end)])
end
xlabel('Time (sec)'); sgtitle('\bfṼ_2','fontsize',14);
set(gcf,'unit','normalized','position',[0.1,0.1,0.4,0.75])
% saveas(gcf,fullfile(SavePath,['Fig13_V2_累積組合.png']))

figure (12)
for j = 1:size(vel_accumu,2)
    subplot(size(vel_accumu,2),1,j)
    plot(Time,cumtrapz(Time,vel_accumu(:,j))*10); grid on
    ylabel(yaxisname{j});xlim([Time(1) Time(end)])
end
xlabel('Time (sec)')
set(gcf,'unit','normalized','position',[0.1,0.1,0.4,0.75])
% saveas(gcf,fullfile(SavePath,['Fig14_V2_累積組合積分.png']))

% [5] D1 + D2
Nt = input("選擇V2累積至IMF_Nt: Nt = "); % 由"EMD累積組合積分"決定
Dt = cumtrapz(Time , vel_accumu(:,Nt))*10;
D2 = D1 + Dt;
figure (13)
plot(Time , D1 , 'LineWidth',1); grid on; hold on ;
plot(Time , D2 , 'LineWidth',1) ; ylabel('Displacement (mm)'); xlabel('Time (sec)');
xlim([Time(1) Time(end)]) ; legend('Dt','D2');
set(gcf,'unit','normalized','position',[0.2,0.2,0.4,0.25]);
% saveas(gcf,fullfile(SavePath,['Fig15_D1&D2(Nt=',num2str(Nt),').png']))
figure;
plot(Time , Dt , 'LineWidth',1) ; ylabel('Displacement (mm)'); xlabel('Time (sec)');
xlim([Time(1) Time(end)]) ; legend('Dt');
set(gcf,'unit','normalized','position',[0.2,0.2,0.4,0.25]);

% v2的Arias intensity及斜率 新增************************************************
AI = cumsum(vel_accumu(:,Nt).^2);
AI = AI./AI(end);
for j = 1:length(AI)-1
    slope(j) = (AI(j+1)-AI(j))*SR;
end
figure (14)
subplot(311)
plot(Time , vel_accumu(:,Nt) , 'LineWidth',1); grid on;
ylabel('Velocity (cm/sec)'); xlabel('Time (sec)'); title('V_2')
xlim([Time(1) Time(end)]) ;

subplot(312)
plot(Time , AI , 'LineWidth',1); grid on;
ylabel('Arias intensity'); xlabel('Time (sec)');
xlim([Time(1) Time(end)]) ;

subplot(313)
plot(Time(1:end-1) , slope , 'LineWidth',1); grid on;
ylabel('Slope'); xlabel('Time (sec)');
xlim([Time(1) Time(end)]) ;
set(gcf,'unit','normalized','position',[0.2,0.2,0.35,0.5]);

if mean(Dt(end-10*SR:end))<0
    [tmp,a] = findpeaks(Dt);
else
    [tmp,a] = findpeaks(-Dt);
end

[~,seat] = max(tmp);

D2 = D2 - (D2(a(seat))-D1(a(seat)));
% tp = Time(a(seat));
tp = 48.29;

figure (15)
plot(Time , D1 , 'LineWidth',1); grid on ; hold on ;
plot(Time , D2 , 'LineWidth',1);
plot(tp , D1(a(seat)) ,'ko', 'MarkerSize', 5 , 'LineWidth',1);
ylabel('Displacement (mm)') ; xlabel('Time(sec)') ;
xlim([Time(1) Time(end)]) ;
legend('D1','Shifted waveform of D_2 (t)',['t_p = ',num2str(tp,'%10.2f'),' sec'],'Location','northeast','fontsize',8);
set(gcf,'unit','normalized','position',[0.2,0.2,0.4,0.25]);
% saveas(gcf,fullfile(SavePath,['Fig16_D1&D2平移.png'])) ; 

% [6] D
D = [D1(1:a(seat));D2(a(seat)+1:end)]; % D1"前面"資料是正確的，D2"後面"資料是正確的
figure (16)
plot(Time , D ,'LineWidth',1); grid on; hold on;
figure (17)
%plot(Time, Reladis, 'LineWidth',1); grid on; hold on;
%plot(Time , LVDT(:,i),'r-','LineWidth',1);
xlim([Time(1) Time(end)]); %title([StationName,' - ',DataSequence{MaxPGArSeat}],'fontsize',14);
ylabel('Displacement (mm)') ; xlabel('Time (sec)'); %ylim([-10 15])
% legend(['Estimated (Max = ',num2str(max(abs(D)),'%10.2f'),' mm / Residual = ',num2str(abs(mean(D(end-5*SR:end))),'%10.2f'),' mm)'],['Measured (Max = ',num2str(max(abs(LVDT(:,i))),'%10.2f'),' mm / Residual = ',num2str(abs(nanmean(LVDT(end-5*SR:end,i))),'%10.2f'),' mm)'],'Location',"southoutside");
set(gcf,'unit','normalized','position',[0.2,0.2,0.4,0.35])
% saveas(gcf,fullfile(SavePath,['Fig17_D_final.png'])) ;



figure (18)
subplot(2,1,1)
plot(Time , Accel); grid on
axis([Time(1) Time(end) max(abs(Accel))*-1.1  max(abs(Accel))*1.1])
ylabel('Acceleration (gal)') ; set(gca,'FontSize',10)
%title([StationName,' - ',DataSequence{MaxPGArSeat}],'fontsize',14);

subplot(2,1,2)
plot(Time , D ,'LineWidth',1); grid on; hold on;
xlim([Time(1) Time(end)]); 
%ylim([-10 20]) ;
% set(gca,'FontSize',10)
ylabel('Displacement (mm)') ; xlabel('Time (sec)'); %ylim([-15 15])
set(gcf,'unit','normalized','position',[0.3,0.3,0.35,0.45]);
% saveas(gcf,fullfile(SavePath,['Fig18_我的結果.png'])) ;

% 與位移計量測結果比較 新增****************************************************
figure(19)
plot(Time,D,'LineWidth',1); grid on; hold on;
%plot(Time,Reladis,'r','LineWidth',1);
xlim([Time(1) Time(end)]);
%ylim([-10 20]) ; 
set(gca,'FontSize',10); legend('Estimated','Measured')
ylabel('Displacement (mm)') ; xlabel('Time (sec)'); %ylim([-15 15])
set(gcf,'unit','normalized','position',[0.3,0.3,0.4,0.35]);



function [yrc] = SSA(SR , x , pt , ch , nrh , nch , ~)
% x : Time Series
% pt : Sample Point(Data Length)
% ch : Number of channel
% nrh : Row number of Hankel Matrix
% nch : column number of Hankel Matrix
% yrc : Reconstruction signal

if size(x,1) <  size(x,2)
    if ~isdeployed
        sprintf('Data check: Data size is Error.')
    end
    return
end


% ------------------
% >>> Step 1 : Embedding
% ------------------
X = zeros(ch*nrh,nch);
for I = 1 : ch
    X(I:ch:end,:) = hankel(x(1:nrh,I),x(nrh:pt,I));
end


% ------------------
% >>> Step 2 : Singular Value Decomposition (SVD)
% ------------------
[U,S,V] = svd(X,'econ');
V = V';
S = diag(S);

figure
plot(S(1:20),'-o') ;  grid on ; legend(['L = ',num2str(nrh)])
xlabel('Number of Singular Value','FontSize',11) ;
ylabel('Singular Value','FontSize',11)
title(['Singular Spectrum'],'FontSize',14)
% set(gcf,'unit','normalized','position',[0.3,0.3,0.35,0.45]);
set(gcf,'unit','normalized','position',[0.2,0.2,0.4,0.25]);
% saveas(gcf,fullfile(SavePath,['Fig3_SV.png'])) ;

% ------------------
% >>> Step 3 : Grouping
% ------------------
E = U(:,2:end) *(diag(S(2:end,1))) * V(2:end,:); % Detail Part : 除了 point 1 以外的其他點


% ------------------
% >>> Step 4 : Reconstruction
% ------------------
yrc = zeros(pt,ch);
for I = 1 : ch
    % Separate Hankel matrix according to each channel
    tmph = E(I:ch:end,:);
    tmph = flipud(tmph);
    % Reconstriction each measurement
    for J = 1 : pt
        tmp = diag(tmph,J-nrh);
        yrc(J,I) = mean(tmp);
    end
end

Time = (0 : 1/SR : (length(yrc)-1)*1/SR)';
AxisValue = max(max(abs(yrc)));

figure
plot(Time,yrc) % Reconstriction signal
axis([Time(1) Time(end) AxisValue*-1.05 AxisValue*1.05])
set(gca,'FontSize',10) ;
ylabel('Accel. (gal)','FontSize',11) ; xlabel('Time (sec)','FontSize',11) ;
title('Reconstructed Acceleration','FontSize',14)
% set(gcf,'unit','normalized','position',[0.3,0.3,0.35,0.45]);
set(gcf,'unit','normalized','position',[0.2,0.2,0.4,0.25]);
% saveas(gcf,fullfile(SavePath,['Fig4_RecAccel.png'])) ;

end


function [Dwmm] = drift_ormsby(Accel , SR , Time , fcut , ft)

% zero mean
Accel = Accel - mean(Accel(1:SR)); %!!!!!!!!!!!!!!!!!!!!!!

% 基線修正
Acc = ormsby(Accel , ft , fcut , SR , length(Accel));
Cor_acc = Acc - mean(Acc(1:2*SR));

vel = cumtrapz(Time, Cor_acc);
vel_2 = detrend(vel,1);
trend = vel - vel_2;

disp_2 = cumtrapz(Time, vel_2);
Dwmm = disp_2*10; %cm-mm   % 位移 mm

end


function [y] = ormsby(x , ft , fcut , SR , DataLenge)
% fcut = 1.8;  % 參數 (從0 - 1.8Hz值為0)
% ft = 2;    % 參數 (1.8 - 2Hz區間進行訊號修正)

dt = 1/SR;
Nw = DataLenge+1;
mnw = fix(Nw/2);
Nww = 1:1:Nw;
w = kaiser(Nw,10);
orm(mnw+1) = 2*(fcut+ft)/2;
t1 = 0;

for nn = mnw + 2:Nw
    t1 = t1 + dt;
    orm(nn) = (cos(2*pi*fcut*t1) - cos(2*pi*ft*t1))/(ft-fcut)/t1^2/2/(pi)^2;
end

for ii = 1:mnw
    orm(ii) = orm(Nw+1-ii);
end

% window time interval
tw = w.*orm';
Xaf = conv(x,tw,'same')*dt;
y = x-Xaf; % minus low pass portion
end

