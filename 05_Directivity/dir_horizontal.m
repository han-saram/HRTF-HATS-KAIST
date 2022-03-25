%-------------------------------------------------------------------------
%   Date : July 22, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : horizontal plane directivity
%   Synopsis : plot horizontal plane directivity
%	Algorithm : -
%-------------------------------------------------------------------------

clc
clear
close all

% elevation angle (0 deg)
elev = 0;

% azimuth angles (-180 ~ +180 deg)
azims = (-180:5:180)';
N_azim = length(azims);

% initialization
Fs = 48e3;
N = 960;
f = (0:N/2-1)*Fs/N;
dir_map_L = zeros(N_azim,N/2);
dir_map_R = zeros(N_azim,N/2);

% reference for directivity
[h_L0,h_R0] = hrir_hats_F(0,elev);
N0 = length(h_R0);
h_L = zeros(N,1);
h_R = zeros(N,1);
h_L(1:N0) = h_L0;
h_R(1:N0) = h_R0;

H_L = abs(fft(h_L));
H_R = abs(fft(h_R));
H_L_ref = H_L(1:N/2);
H_R_ref = H_R(1:N/2);

% directivity map
for i = 1:N_azim
    azim = azims(i);
    
    [h_L0,h_R0] = hrir_hats_F(azim,elev);
    h_L = zeros(N,1);
    h_R = zeros(N,1);
    h_L(1:N0) = h_L0;
    h_R(1:N0) = h_R0;
    
    H_L = abs(fft(h_L));
    H_R = abs(fft(h_R));
    H_L = H_L(1:N/2);
    H_R = H_R(1:N/2);
    
    dir_map_L(i,:) = 20*log10(H_L./H_L_ref);
    dir_map_R(i,:) = 20*log10(H_R./H_R_ref);
end

% plot
max = 10;
min = -30;
tic = ceil((max-min)/5);
lims = [max min tic];

figure
dirplot_F(azims,dir_map_L(:,f == 0.75e3),'-r',lims);
hold on
dirplot_F(azims,dir_map_L(:,f == 1.5e3),'-g');
hold on
dirplot_F(azims,dir_map_L(:,f == 3.0e3),'-c');
hold on
dirplot_F(azims,dir_map_L(:,f == 6.0e3),'-b');
hold on
dirplot_F(azims,dir_map_L(:,f == 12.0e3),'-m');
legend('0.75 kHz','1.5 kHz','3.0 kHz','6.0 kHz','12.0 kHz');
title('Left ear horizontal plane directivity','fontsize',18);
set(gcf,'position',[100 500 800 800]);

figure
dirplot_F(azims,dir_map_R(:,f == 0.75e3),'-r',lims);
hold on
dirplot_F(azims,dir_map_R(:,f == 1.5e3),'-g');
hold on
dirplot_F(azims,dir_map_R(:,f == 3.0e3),'-c');
hold on
dirplot_F(azims,dir_map_R(:,f == 6.0e3),'-b');
hold on
dirplot_F(azims,dir_map_R(:,f == 12.0e3),'-m');
legend('0.75 kHz','1.5 kHz','3.0 kHz','6.0 kHz','12.0 kHz');
title('Right ear horizontal plane directivity','fontsize',18);
set(gcf,'position',[1000 500 800 800]);
