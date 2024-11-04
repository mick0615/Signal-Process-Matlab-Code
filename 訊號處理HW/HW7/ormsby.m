%% ormsby function
% Low-pass Ormsby filter

% orm(f) = 2[cos(2*pi*fc*t)-cos(2*pi*ft*t)]/(ft-fc)/t^2
% fcut and ft are corner and transient frequency

function [y] = ormsby(x,ft,fcut,pt,fs,d_t)
%% Create the Kaiser window
% low pass filter after baseline correlation
% apply Ormsby low-pass filter with Kaiser window
% call window function
% -------------------------------------------------------------------------
Nw = pt;
mnw = fix(Nw/2);
Nww = 1:1:Nw;
w = kaiser(Nw,10); % w = window(Nw, 'kais', 10) ; % kaiser window
orm(mnw+1) = 2*(fcut+ft)/2; % at t=0 amplitude is 2fa = 2*(fc+ft)/2
t1=0;

for nn=mnw+2:Nw
    t1=t1+d_t;
    orm(nn)=(cos(2*pi*fcut*t1)-cos(2*pi*ft*t1))/(ft-fcut)/t1^2/2/(pi)^2;
end
for ii=1:mnw
    orm(ii)=orm(Nw+1-ii);
end
%window time interval
tw=w.*orm'; 
Xaf=conv(x,tw,'same')*d_t; % y = Xaf (lowpass)
y = x - Xaf;              % y = x-Xaf (high pass)
end