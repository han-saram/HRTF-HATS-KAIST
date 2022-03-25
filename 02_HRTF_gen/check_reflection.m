%-------------------------------------------------------------------------
%   Date : June 30, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------

clc
clear
close all

% elevation angles (-40 ~ +90 deg)
elevs = (-40:5:90)';
N_elev = length(elevs);

% plot limit
h_lims = [-0.05 0.05];

%% azimuth X deg, elevation 0 deg (horizontal plane)
% source direction
azim = 90;
elev = 0;

% root of directory
root = '../00_Data/TF_HATS';

% load IR
[h,Fs] = ir_hats_F(azim,elev,root);

% time axis
N = length(h);
time_axis = 1000*(0:N-1)/Fs;            % ms

% plot
figure
plot(time_axis,h(:,1));
hold on
plot(time_axis,h(:,2));
legend('Left Ear','Right Ear');
xlabel('Time [ms]','fontsize',12); ylabel('Amplitude','fontsize',12);
title(['IR amplitude ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([0 10 h_lims(1) h_lims(2)]);
grid on
set(gcf,'position',[100 500 800 400]);

%% azimuth X deg, elevation -40 ~ +90 deg
IR_L = zeros(N_elev,N);
IR_R = zeros(N_elev,N);

% elevation -40 ~ +90 deg
for elev_idx = 1:N_elev
    elev = elevs(elev_idx);
    
    % load IR
    [h,~] = ir_hats_F(azim,elev,root);
    IR_L(elev_idx,:) = h(:,1);
    IR_R(elev_idx,:) = h(:,2);
end

% plot
figure
imagesc(time_axis,elevs,IR_L,h_lims); axis xy
colormap jet
xlabel('Time [ms]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
title(['Left ear IR amplitude ( Azimuth ',num2str(azim),' deg )'],'fontsize',14);
axis([0 10 elevs(1) elevs(end)]);
grid on
set(gcf,'position',[900 500 800 400]);

figure
imagesc(time_axis,elevs,IR_R,h_lims); axis xy
colormap jet
xlabel('Time [ms]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
title(['Right ear IR amplitude ( Azimuth ',num2str(azim),' deg )'],'fontsize',14);
axis([0 10 elevs(1) elevs(end)]);
grid on
set(gcf,'position',[1700 500 800 400]);
