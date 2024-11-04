% M11205314-HW2-2-Code
clc; clear; close all;
set(0, 'DefaultAxesFontName', 'times');
set(0,'DefaultAxesFontSize',10);
data = xlsread('data.xlsx');
timeAxis = data(:,1);
data5 = data(:,3);
data5 = data5 - mean(data5);
dt = timeAxis(2) - timeAxis(1);
fs = round(1/dt); % sampling rate
Ndata = length(timeAxis); % number of points

% window function
wsigmaList = [0.4, 1.0, 4.0]; % sigma for the Gausian Window
wlen = 5; % sec, window length
wtime = (-wlen:dt:wlen)'; % window time
wdata = zeros(length(wtime),1); % data in the window
wNdata = length(wtime); % window number of points
freqAxis = fs*((1:wNdata)'-1)/(wNdata-1);

phySpec = zeros(wNdata, Ndata,length(wsigmaList)); % physical spectrum
for isig = 1:length(wsigmaList)
    wsigma = wsigmaList(isig);
    wfun = sqrt(2/(sqrt(2*pi)*wsigma))*exp(-wtime.^2/wsigma^2); % Gausian window function
    for I = 1:Ndata
        idx1 = I-(wNdata-1)/2;
        idx2 = idx1 + (wNdata-1);
        if idx1 <= 0
            % symmetry on the edge
            wdata = [flipud(data5(1:(wNdata-idx2))); data5(1:idx2)];
            wdata = wfun.*wdata;
        elseif idx1 > (Ndata-wNdata+1)
            % symmetry on the edge
            Nflip = wNdata - (Ndata - idx1 + 1);
            wdata = [data5(idx1:end); flipud(data5((Ndata-Nflip+1):end))];
            wdata = wfun.*wdata;
        else
            % normal procedure
            wdata = wfun.*data5(idx1:idx2);
        end
        phySpec(:,I,isig) = abs(fft(wdata))/wNdata*2;
    end
    wfunStore(:,isig) = wfun;
end
%% plot the result(time response)
plot(timeAxis,data5); title('Response Data Signal'); xlabel('Time (sec)'); ylabel('Response (g)'); grid on; grid minor;
set(gcf,'Position',[50,50,600,200]);
%saveas(gcf, folderLocation+"\"+plotFolder+"\Response Data", 'png')
%saveas(gcf, folderLocation+"\"+plotFolder+"\Response Data", 'fig') 
%% plot the result(window function)
figure()
plot(wtime,wfunStore)
set(gcf,'Position',[50,50,400,200])
grid on; grid minor
xlabel("Time (sec)")
legend("\sigma = 0.4", "\sigma = 1.0", "\sigma = 4.0")
%% plot the result(physical spectrum)
for isig = 1:length(wsigmaList)
    figure()
    set(gcf,'Position',[50,50,700,500])
    surf(timeAxis,freqAxis,phySpec(:,:,isig),'LineStyle','none','FaceLighting','phong');
    title("Synthetic Time Series Physical Spectrum"+", Gausian Window with \sigma = "+wsigmaList(isig));
    colormap('jet')
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    xlim([0,40])
    ylim([0,20])
    view(0,90)
    colorbar
end
%     saveas(gcf, folderLocation+"\"+plotFolder+"\Synthetic Time series Physical Spectrum sigma = "+isig, 'png')
%     saveas(gcf, folderLocation+"\"+plotFolder+"\Synthetic Time series Physical Spectrum sigma = "+isig, 'fig')
