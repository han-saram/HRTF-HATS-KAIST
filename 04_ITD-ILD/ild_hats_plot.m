%-------------------------------------------------------------------------
%   Date : July 16, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : HATS ILD profile & map
%   Synopsis : plot HATS ILD profile & map
%	Algorithm : -
%-------------------------------------------------------------------------

clc
clear
close all

% elevation angles (-40 ~ +90 deg)
elev = 0;

% azimuth angles (-180 ~ +180 deg)
azims = (-180:5:180)';
N_azim = length(azims);

% initialization
N_f = 2400;
ILD_map = zeros(N_f,N_azim);

% ILD Map
for j = 1:N_azim
    azim = azims(j);
    
    % ILD retrieval
    [ILD,f] = ild_hats_F(azim,elev);
    ILD_map(:,j) = ILD;
end

% plot
figure
plot(azims,ILD_map(f == 800,:),'LineWidth',1.2);
hold on
plot(azims,ILD_map(f == 1600,:),'LineWidth',1.2);
hold on
plot(azims,ILD_map(f == 3200,:),'LineWidth',1.2);
hold on
plot(azims,ILD_map(f == 6400,:),'LineWidth',1.2);
hold on
plot(azims,ILD_map(f == 12800,:),'LineWidth',1.2);
legend('0.8 kHz','1.6 kHz','3.2 kHz','6.4 kHz','12.8 kHz','Location','northwest');
xlabel('Azimuth [deg.]','fontsize',12); ylabel('Level [dB]','fontsize',12);
title('ILD profile of HATS in the horizontal plane','fontsize',14);
axis([azims(1) azims(end) -40 40]);
grid on
set(gcf,'position',[100 500 700 500]);

figure
imagesc(azims,0.001*f,ILD_map); axis xy
colormap jet; colorbar
xlabel('Azimuth [deg.]','fontsize',12); ylabel('Frequency [kHz]','fontsize',12);
title('ILD map of HATS in the horizontal plane','fontsize',14);
title(colorbar,'Level [dB]')
axis([azims(1) azims(end) 0.001*f(1) 0.001*f(end)]);
grid on
set(gcf,'position',[900 500 700 500]);
