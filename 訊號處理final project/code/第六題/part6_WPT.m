clc; clear; close all;
%% Data
Data = readmatrix("Ming-Li School 1999-9-21 data.xlsx");

% For analysis define the duration starting from 20 sec to the end record.
F8 = Data(4001:end,8);
F15 = Data(4001:end,15);
F21 = Data(4001:end,21);
F28 = Data(4001:end,28);

SR = 200;
dt = 1/SR;
N = length(F28);
t = 0:dt:(N-1)*dt;

%% Low-Pass and Zero-Mean
cff = 5;
[b,a] = butter(5,cff/(0.5*SR),'low');

F21 = filtfilt(b,a,F21);
F21 = F21-mean(F21);
F8 = filtfilt(b,a,F8);
F8 = F8-mean(F8);
F15 = filtfilt(b,a,F15);
F15 = F15-mean(F15);
F28 = filtfilt(b,a,F28);
F28 = F28-mean(F28);

%% Bior6.8,Level-9
u = F21;
Level = 9;

%     [tvec,Fvec,Eenergy,Amp,WPTree,RconstCoef,Fint] = fn_wptspect(u,SR,Level);

if size(u,1)==1
    u=u';
end
wname = 'bior6.8'; % selected wavelet base
dt = 1/SR; % dt
N = max(size(u)); % number of data
tvec = [0:dt:(N-1)*dt]'; % time

% Index of Paley order
index = 0;
for jlevel=1:Level
    index = [index fliplr(index)+2^(jlevel-1)]; % sift index (component index)
end
index = index+1;

% Wavelet Packet Transform: wpdec()=WPT(Wavelet Packet Decomposition)
WPTree = wpdec(u,Level,wname); % wavelet packet tree
% plot(WPTree)

Amp = zeros(length(u),2^Level);
for icomp = 1 : 1 : 2^Level
    % Reconstruction of all components
    RconstCoef(:,icomp)= wprcoef(WPTree,[Level,index(icomp)-1]); % wprcoef()=Wavelet Packet "Reconstruct" Coefficient

    % Hibert amplitude of all components (for calculating Energy)
    Amp(:,icomp) = abs(hilbert(RconstCoef(:,icomp)));
end
% Energy = square of Hilbert amplitude
Energy = Amp.*Amp;

E1 = Energy;
% Frequency resolution
Fvec1 = linspace(0,1/dt/2,2^Level)';

% Interested frequency
Fint1 = find(Fvec1<50);

% plot time-freq spectrogram
figure
surf(tvec,Fvec1(Fint1),E1(:,Fint1)','LineStyle','none','FaceColor','interp','FaceLighting','phong','EdgeColor','interp')
view(0,90),colormap('jet'),ylim([0 5]),xlabel('Time (sec)'),ylabel('Frequency (Hz)'),xlim([0 tvec(end)]);
colorbar
title("Bior6.8,Level-9")

%% SSA
F = [F21];

[pt,ch]=size(F);
totalpair = 7;
nrh = 300;
nch = pt-nrh+1;
X = zeros(ch*nrh,nch);

for i = 1:ch
    X(i:ch:end,:) = hankel(F(1:nrh,i),F(nrh:pt,i));
end

if true
 [U,S,V] = svd(X,'econ');
else
 Scov = X*X';
 [U,S,V] = svd(Scov,'econ');
 clear Scov
end

S = diag(S);

figure
plot(1:length(S),S,'-o','LineWidth',1.5)
xlim([0,20])
xlabel('Singular Value')
ylabel('Singular Value')
title('Singular Value Plot')

%%
X_recon = zeros([size(X),totalpair]);

for i = 1:totalpair
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

x_recon = zeros(pt,ch,totalpair);
for K = 1:totalpair
    for I = 1:ch
        tmph = X_recon(I:ch:end,:,K);
        tmph = flipud(tmph);
        for J= 1:pt
            tmp=diag(tmph,J-nrh);
            x_recon(J,I,K) = mean(tmp);
        end
    end
end
clear tmph

ReCon_1 = x_recon(:,1,1) + x_recon(:,1,2) + x_recon(:,1,3) + x_recon(:,1,4) + x_recon(:,1,5) + x_recon(:,1,6) + x_recon(:,1,7);

figure;
subplot(7,1,1)
plot(t,F(:,1),'LineWidth',1.5)
hold on
plot(t,ReCon_1,'LineWidth',1.5)
xlabel('Time(s)');
ylabel('Acc.');
legend('Original','Reconstructed');
for i=1:totalpair
    subplot(8,1,i+1)
    plot(t,x_recon(:,1,i),'LineWidth',1.5)
    ylabel('Acc.');
    legend(['Group',sprintf('%d',i)]);
end
sgtitle('CH21 Response');

figure
for i = 1:7
    xi = x_recon(:,1,i);
    Fi = F(:,1);
    XT = fft(xi,N);
    FT = fft(Fi,N);
    magX = abs(XT);
    magF = abs(FT);
    f = (0:N-1)*SR/N;
    subplot(7,1,i)
    hold on
    plot(f(1:N/2),magX(1:N/2)*2/N,'b','LineWidth',1.5);
    plot(f(1:N/2),magF(1:N/2)*2/N,'r--','LineWidth',1.5);
    legend('reconstructed','original')
    xlabel('Frequency(Hz)')
    ylabel('amplitude')
    xlim([0 5])
    title(['Group',num2str(i)])
    hold off
end

%% Bior6.8,Level-9
u = x_recon(:,1,1);
Level = 9;

%     [tvec,Fvec,Eenergy,Amp,WPTree,RconstCoef,Fint] = fn_wptspect(u,SR,Level);

if size(u,1)==1
    u=u';
end
wname = 'bior6.8'; % selected wavelet base
dt = 1/SR; % dt
N = max(size(u)); % number of data
tvec = [0:dt:(N-1)*dt]'; % time

% Index of Paley order
index = 0;
for jlevel=1:Level
    index = [index fliplr(index)+2^(jlevel-1)]; % sift index (component index)
end
index = index+1;

% Wavelet Packet Transform: wpdec()=WPT(Wavelet Packet Decomposition)
WPTree = wpdec(u,Level,wname); % wavelet packet tree
% plot(WPTree)

Amp = zeros(length(u),2^Level);
for icomp = 1 : 1 : 2^Level
    % Reconstruction of all components
    RconstCoef(:,icomp)= wprcoef(WPTree,[Level,index(icomp)-1]); % wprcoef()=Wavelet Packet "Reconstruct" Coefficient

    % Hibert amplitude of all components (for calculating Energy)
    Amp(:,icomp) = abs(hilbert(RconstCoef(:,icomp)));
end
% Energy = square of Hilbert amplitude
Energy = Amp.*Amp;

E2 = Energy;
% Frequency resolution
Fvec2 = linspace(0,1/dt/2,2^Level)';

% Interested frequency
Fint2 = find(Fvec2<50);

% plot time-freq spectrogram
figure
surf(tvec,Fvec2(Fint2),E2(:,Fint2)','LineStyle','none','FaceColor','interp','FaceLighting','phong','EdgeColor','interp')
view(0,90),colormap('jet'),ylim([0 5]),xlabel('Time (sec)'),ylabel('Frequency (Hz)'),xlim([0 tvec(end)]);
colorbar
title("Bior6.8,Level-9")

%%
figure
subplot(211)
surf(tvec,Fvec1(Fint1),E1(:,Fint1)','LineStyle','none','FaceColor','interp','FaceLighting','phong','EdgeColor','interp')
view(90,0),colormap('jet'),ylim([0 5]),zlim([0 12000]),xlabel('Time (sec)'),ylabel('Frequency (Hz)');zlabel('Energy'),xlim([0 tvec(end)]);
title('Roof Response (distribution of component energy)')
colorbar
subplot(212)
surf(tvec,Fvec2(Fint2),E2(:,Fint2)','LineStyle','none','FaceColor','interp','FaceLighting','phong','EdgeColor','interp')
view(90,0),colormap('jet'),ylim([0 5]),zlim([0 12000]),xlabel('Time (sec)'),ylabel('Frequency (Hz)');zlabel('Energy'),xlim([0 tvec(end)]);
title('Roof Response(first mode modal contribution)')
colorbar
%% Calculate First Mode Contribution
Contri = sum(sum(E2(:,Fint2)))/sum(sum(E1(:,Fint1)));
Contri = Contri*100;
disp(['First Mode Contribution = ',num2str(Contri),'%'])