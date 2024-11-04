clc;clear;close all;
data = xlsread('HW8-data');
t4_base_acc = data(:,2);
t4_top_acc = data(:,3);
t5_base_acc = data(:,4);
t5_top_acc = data(:,5);
t7_base_acc = data(:,6);
t7_top_acc = data(:,7);

t4_acc = t4_top_acc-t4_base_acc;
t5_acc = t5_top_acc-t5_base_acc;
t7_acc = t7_top_acc-t7_base_acc;


SR = 256;
delt = 1/SR;
N = 10241;
t = 0:delt:(N-1)*delt;
duration = N*delt;
ndata = duration*SR+1;

% fcut = 0.1;
% ft = 0.15;
% t4_acc = ormsby(t4_acc,ft,fcut,length(t),SR,delt);
% t5_acc = ormsby(t5_acc,ft,fcut,length(t),SR,delt);
% t7_acc = ormsby(t7_acc,ft,fcut,length(t),SR,delt);


t4_acc = t4_acc - mean(t4_acc);
t5_acc = t5_acc - mean(t5_acc);
t7_acc = t7_acc - mean(t7_acc);



hb_t4_acc = hilbert(t4_acc);
hb_t5_acc = hilbert(t5_acc);
hb_t7_acc = hilbert(t7_acc);

%% Method 1(test4)
tmpr4 = real(hb_t4_acc);
tmpi4 = imag(hb_t4_acc);
    IF_Ufh1_4 = (tmpr4.*([0;diff(tmpi4)])-tmpi4.*([0;diff(tmpr4)]))./(tmpr4.^2+tmpi4.^2);
    IF_Ufh14(:) = IF_Ufh1_4*SR/2/pi;
figure;
hold on
plot(t,IF_Ufh14);ylim([-400 400]);xlabel('t(sec)');ylabel('IF(Hz)');title('test4 (Method 1)')
%% Method 2(test4)
IF_Ufh24 = zeros(N-1,1);
Ufh24 = hb_t4_acc;
IF_Ufh24 = diff(unwrap(angle(Ufh24)))*SR/2/pi;
t1 = 0:delt:(N-2)*delt;
plot(t1,IF_Ufh24);xlabel('t(sec)');ylabel('IF(Hz)');title('test4 (Method 2)')
%% Method 3(test4)
IF_t34 = zeros(ndata,1);
    Ufh34 = hb_t4_acc;
    theta = unwrap(angle(Ufh34));
    for idata = 1:(ndata-1)
        if or(idata == 1, idata == 2)
            IF_Ufh34(idata) = (-25*theta(idata)+48*theta(idata+1)-36*theta(idata+2)+16*theta(idata+3)-3*theta(idata+4))/12;
        elseif or(idata == ndata-1, idata == ndata-2)
            IF_Ufh34(idata) = (+25*theta(idata)-48*theta(idata-1)+36*theta(idata-2)-16*theta(idata-3)+3*theta(idata-4))/12;
        else
            IF_Ufh34(idata) = (  theta(idata-2) -8*theta(idata-1) +8*theta(idata+1)   -theta(idata+2)               )/12;
        end
    end
    IF_Ufh34(:) = IF_Ufh34(:)*SR/2/pi;
plot(t,IF_Ufh34);xlabel('t(sec)');ylabel('IF(Hz)');title('Instantaneous Frequency,IF(test4)')
legend('Method 1:Direct Differentiation','Method 2:Derivative of Arctan','Method 3:Finite Differentiation')
hold off
%% Method 4(test4)
figure;
IF_Ufh44(:)=instfreq(t4_acc,SR);
instfreq(t4_acc,SR);title('Instantaneous Frequency,IF(test4) Method 4:instfreq')
%% smooth all of them
Naverage = 9;

    IF_Ufh14_smooth(:) = movmean(IF_Ufh14(:),Naverage);
    IF_Ufh24_smooth(:) = movmean(IF_Ufh24(:),Naverage);
    IF_Ufh34_smooth(:) = movmean(IF_Ufh34(:),Naverage);
    IF_Ufh44_smooth(:) = movmean(IF_Ufh44(:),Naverage);

for i = 1:100

        IF_Ufh14_smooth(:) = movmean(IF_Ufh14_smooth(:),Naverage);
        IF_Ufh24_smooth(:) = movmean(IF_Ufh24_smooth(:),Naverage);
        IF_Ufh34_smooth(:) = movmean(IF_Ufh34_smooth(:),Naverage);
        IF_Ufh44_smooth(:) = movmean(IF_Ufh44_smooth(:),Naverage);
end
figure;
hold on
plot(t,detrend(IF_Ufh14_smooth,10));
plot(t1,detrend(IF_Ufh24_smooth,10));
plot(t,detrend(IF_Ufh34_smooth,10));xlabel('t(sec)');ylabel('IF(Hz)');title('Instantaneous Frequency(smooth),IF(test4)')
legend('Method 1:Direct Differentiation','Method 2:Derivative of Arctan','Method 3:Finite Differentiation')
hold on
%% test 4 (subplot)
figure;
subplot(3,1,1);plot(t,IF_Ufh14);xlabel('t(sec)');ylabel('IF(Hz)');title('test4 (Method 1)');
subplot(3,1,2);plot(t1,IF_Ufh24);xlabel('t(sec)');ylabel('IF(Hz)');title('test4 (Method 2)');
subplot(3,1,3);plot(t,IF_Ufh34);xlabel('t(sec)');ylabel('IF(Hz)');title('test4 (Method 3)');
figure;
subplot(3,1,1);plot(t,detrend(IF_Ufh14_smooth,10));xlabel('t(sec)');ylabel('IF(Hz)');title('test4 smooth (Method 1)')
subplot(3,1,2);plot(t1,detrend(IF_Ufh24_smooth,10));xlabel('t(sec)');ylabel('IF(Hz)');title('test4 smooth (Method 2)')
subplot(3,1,3);plot(t,detrend(IF_Ufh34_smooth,10));xlabel('t(sec)');ylabel('IF(Hz)');title('test4 smooth (Method 3)')
%% Method 1(test5)
tmpr5 = real(hb_t5_acc);
tmpi5 = imag(hb_t5_acc);
    IF_Ufh1_5 = (tmpr5.*([0;diff(tmpi5)])-tmpi5.*([0;diff(tmpr5)]))./(tmpr5.^2+tmpi5.^2);
    IF_Ufh15(:) = IF_Ufh1_5*SR/2/pi;
figure;
hold on
plot(t,IF_Ufh15);ylim([-400 400]);xlabel('t(sec)');ylabel('IF(Hz)');title('test5 (Method 1)')
%% Method 2(test5)
IF_Ufh25 = zeros(N-1,1);
Ufh25 = hb_t5_acc;
IF_Ufh25 = diff(unwrap(angle(Ufh25)))*SR/2/pi;
t1 = 0:delt:(N-2)*delt;
plot(t1,IF_Ufh25);xlabel('t(sec)');ylabel('IF(Hz)');title('test5 (Method 2)')
%% Method 3(test5)
IF_t35 = zeros(ndata,1);
    Ufh35 = hb_t5_acc;
    theta = unwrap(angle(Ufh35));
    for idata = 1:(ndata-1)
        if or(idata == 1, idata == 2)
            IF_Ufh35(idata) = (-25*theta(idata)+48*theta(idata+1)-36*theta(idata+2)+16*theta(idata+3)-3*theta(idata+4))/12;
        elseif or(idata == ndata-1, idata == ndata-2)
            IF_Ufh35(idata) = (+25*theta(idata)-48*theta(idata-1)+36*theta(idata-2)-16*theta(idata-3)+3*theta(idata-4))/12;
        else
            IF_Ufh35(idata) = (  theta(idata-2) -8*theta(idata-1) +8*theta(idata+1)   -theta(idata+2)               )/12;
        end
    end
    IF_Ufh35(:) = IF_Ufh35(:)*SR/2/pi;
plot(t,IF_Ufh35);xlabel('t(sec)');ylabel('IF(Hz)');title('Instantaneous Frequency,IF(test5)')
legend('Method 1:Direct Differentiation','Method 2:Derivative of Arctan','Method 3:Finite Differentiation')
hold off
%% Method 4(test5)
figure;
IF_Ufh45(:)=instfreq(t5_acc,SR);
instfreq(t5_acc,SR);title('Instantaneous Frequency,IF(test5) Method 4:instfreq')
%% smooth all of them
Naverage = 9;

    IF_Ufh15_smooth(:) = movmean(IF_Ufh15(:),Naverage);
    IF_Ufh25_smooth(:) = movmean(IF_Ufh25(:),Naverage);
    IF_Ufh35_smooth(:) = movmean(IF_Ufh35(:),Naverage);
    IF_Ufh45_smooth(:) = movmean(IF_Ufh45(:),Naverage);

for i = 1:100

        IF_Ufh15_smooth(:) = movmean(IF_Ufh15_smooth(:),Naverage);
        IF_Ufh25_smooth(:) = movmean(IF_Ufh25_smooth(:),Naverage);
        IF_Ufh35_smooth(:) = movmean(IF_Ufh35_smooth(:),Naverage);
        IF_Ufh45_smooth(:) = movmean(IF_Ufh45_smooth(:),Naverage);
end
figure;
hold on
plot(t,detrend(IF_Ufh15_smooth,10));
plot(t1,detrend(IF_Ufh25_smooth,10));
plot(t,detrend(IF_Ufh35_smooth,10));;xlabel('t(sec)');ylabel('IF(Hz)');title('Instantaneous Frequency(smooth),IF(test5)')
legend('Method 1:Direct Differentiation','Method 2:Derivative of Arctan','Method 3:Finite Differentiation')
hold on
%% test5 subplot
figure;
subplot(3,1,1);plot(t,IF_Ufh15);xlabel('t(sec)');ylabel('IF(Hz)');title('test5 (Method 1)');
subplot(3,1,2);plot(t1,IF_Ufh25);xlabel('t(sec)');ylabel('IF(Hz)');title('test5 (Method 2)');
subplot(3,1,3);plot(t,IF_Ufh35);xlabel('t(sec)');ylabel('IF(Hz)');title('test5 (Method 3)');
figure;
subplot(3,1,1);plot(t,detrend(IF_Ufh15_smooth,10));xlabel('t(sec)');ylabel('IF(Hz)');title('test5 smooth (Method 1)')
subplot(3,1,2);plot(t1,detrend(IF_Ufh25_smooth,10));xlabel('t(sec)');ylabel('IF(Hz)');title('test5 smooth (Method 2)')
subplot(3,1,3);plot(t,detrend(IF_Ufh35_smooth,10));xlabel('t(sec)');ylabel('IF(Hz)');title('test5 smooth (Method 3)')
%% Method 1(test7)
tmpr7 = real(hb_t7_acc);
tmpi7 = imag(hb_t7_acc);
    IF_Ufh1_7 = (tmpr7.*([0;diff(tmpi7)])-tmpi7.*([0;diff(tmpr7)]))./(tmpr7.^2+tmpi7.^2);
    IF_Ufh17(:) = IF_Ufh1_7*SR/2/pi;
figure;
hold on
plot(t,IF_Ufh17);ylim([-400 400]);xlabel('t(sec)');ylabel('IF(Hz)');title('test7 (Method 1)')
%% Method 2(test7)
IF_Ufh27 = zeros(N-1,1);
Ufh27 = hb_t7_acc;
IF_Ufh27 = diff(unwrap(angle(Ufh27)))*SR/2/pi;
t1 = 0:delt:(N-2)*delt;
plot(t1,IF_Ufh27);xlabel('t(sec)');ylabel('IF(Hz)');title('test7 (Method 2)')
%% Method 3(test7)
IF_t37 = zeros(ndata,1);
    Ufh37 = hb_t7_acc;
    theta = unwrap(angle(Ufh37));
    for idata = 1:(ndata-1)
        if or(idata == 1, idata == 2)
            IF_Ufh37(idata) = (-25*theta(idata)+48*theta(idata+1)-36*theta(idata+2)+16*theta(idata+3)-3*theta(idata+4))/12;
        elseif or(idata == ndata-1, idata == ndata-2)
            IF_Ufh37(idata) = (+25*theta(idata)-48*theta(idata-1)+36*theta(idata-2)-16*theta(idata-3)+3*theta(idata-4))/12;
        else
            IF_Ufh37(idata) = (  theta(idata-2) -8*theta(idata-1) +8*theta(idata+1)   -theta(idata+2)               )/12;
        end
    end
    IF_Ufh37(:) = IF_Ufh37(:)*SR/2/pi;
plot(t,IF_Ufh37);xlabel('t(sec)');ylabel('IF(Hz)');title('Instantaneous Frequency,IF(test7)')
legend('Method 1:Direct Differentiation','Method 2:Derivative of Arctan','Method 3:Finite Differentiation')
hold off
%% Method 4(test7)
figure;
IF_Ufh47(:)=instfreq(t7_acc,SR);
instfreq(t7_acc,SR);title('Instantaneous Frequency,IF(test7) Method 4:instfreq')
%% smooth all of them
Naverage = 9;

    IF_Ufh17_smooth(:) = movmean(IF_Ufh17(:),Naverage);
    IF_Ufh27_smooth(:) = movmean(IF_Ufh27(:),Naverage);
    IF_Ufh37_smooth(:) = movmean(IF_Ufh37(:),Naverage);
    IF_Ufh47_smooth(:) = movmean(IF_Ufh47(:),Naverage);

for i = 1:100

        IF_Ufh17_smooth(:) = movmean(IF_Ufh17_smooth(:),Naverage);
        IF_Ufh27_smooth(:) = movmean(IF_Ufh27_smooth(:),Naverage);
        IF_Ufh37_smooth(:) = movmean(IF_Ufh37_smooth(:),Naverage);
        IF_Ufh47_smooth(:) = movmean(IF_Ufh47_smooth(:),Naverage);
end
figure;
hold on
plot(t,detrend(IF_Ufh17_smooth,10));
plot(t1,detrend(IF_Ufh27_smooth,10));
plot(t,detrend(IF_Ufh37_smooth,10));;xlabel('t(sec)');ylabel('IF(Hz)');title('Instantaneous Frequency(smooth),IF(test7)')
legend('Method 1:Direct Differentiation','Method 2:Derivative of Arctan','Method 3:Finite Differentiation')
hold on
%% test 7 subplot
figure;
subplot(3,1,1);plot(t,IF_Ufh17);xlabel('t(sec)');ylabel('IF(Hz)');title('test7 (Method 1)');
subplot(3,1,2);plot(t1,IF_Ufh27);xlabel('t(sec)');ylabel('IF(Hz)');title('test7 (Method 2)');
subplot(3,1,3);plot(t,IF_Ufh37);xlabel('t(sec)');ylabel('IF(Hz)');title('test7 (Method 3)');
figure;
subplot(3,1,1);plot(t,detrend(IF_Ufh17_smooth,10));xlabel('t(sec)');ylabel('IF(Hz)');title('test7 smooth (Method 1)')
subplot(3,1,2);plot(t1,detrend(IF_Ufh27_smooth,10));xlabel('t(sec)');ylabel('IF(Hz)');title('test7 smooth (Method 2)')
subplot(3,1,3);plot(t,detrend(IF_Ufh37_smooth,10));xlabel('t(sec)');ylabel('IF(Hz)');title('test7 smooth (Method 3)')
%% Calculate unwrapped phase
p4 = unwrap(angle(hb_t4_acc));
p5 = unwrap(angle(hb_t5_acc));
p7 = unwrap(angle(hb_t7_acc));
unwraped_phase_test4 = sum(p4)
unwraped_phase_test5 = sum(p5)
unwraped_phase_test7 = sum(p7)
