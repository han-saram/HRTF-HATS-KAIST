%-------------------------------------------------------------------------
%   Date : June 09, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : HATS TF according to elevation
%   Synopsis : plot time & freq. responses of HATS TF measurements
%	Algorithm : -
%-------------------------------------------------------------------------

clc
clear
close all

% azimuth angles (-180 ~ +180 deg)
azims = (-180:5:180)';
N_azim = length(azims);

% phase flag (0: wrap, 1: unwrap)
pha_flag = 1;

% plot limit
h_lims = [-0.05 0.05];          % amplitude
m_lims = [40 100];              % magnitude [dB]
p_lims = [-1 1];                % phase/pi [rad/rad]

%% azimuth 0 deg, elevation -X & X deg (median plane)
% source direction
azim = 0;
elev = 10;              % 0 ~ 40 deg

% root of directory
root = '../00_Data/TF_HATS';

% load IR
U_elev = elev;
D_elev = -elev;
[U_h,Fs] = ir_hats_F(azim,U_elev,root);
[D_h,~] = ir_hats_F(azim,D_elev,root);

% length of IR
N = length(U_h);

% time & freq. axis
time_axis = 1000*(0:N-1)/Fs;            % ms
freq_axis = 0.001*(0:N/2-1)*Fs/N;       % kHz

% FFT
U_H = fft(U_h);
D_H = fft(D_h);

% magnitude
U_TF_mag_L = 20*log10(abs(U_H(1:N/2,1))/20e-6);
U_TF_mag_R = 20*log10(abs(U_H(1:N/2,2))/20e-6);
D_TF_mag_L = 20*log10(abs(D_H(1:N/2,1))/20e-6);
D_TF_mag_R = 20*log10(abs(D_H(1:N/2,2))/20e-6);

% phase
if pha_flag == 0                        % rad/rad
    U_TF_pha_L = angle(U_H(1:N/2,1))/pi;
    U_TF_pha_R = angle(U_H(1:N/2,2))/pi;
    D_TF_pha_L = angle(D_H(1:N/2,1))/pi;
    D_TF_pha_R = angle(D_H(1:N/2,2))/pi;
else
    U_TF_pha_L = unwrap(angle(U_H(1:N/2,1)))/pi;
    U_TF_pha_R = unwrap(angle(U_H(1:N/2,2)))/pi;
    D_TF_pha_L = unwrap(angle(D_H(1:N/2,1)))/pi;
    D_TF_pha_R = unwrap(angle(D_H(1:N/2,2)))/pi;
    p_lims = [min([min(U_TF_pha_L), min(U_TF_pha_R), min(D_TF_pha_L), min(D_TF_pha_R)]), 0];
end

% plot
figure
subplot(321)
plot(time_axis,D_h(:,1));
hold on
plot(time_axis,D_h(:,2));
legend('Left Ear','Right Ear');
xlabel('Time [ms]','fontsize',12); ylabel('Amplitude','fontsize',12);
title(['IR amplitude ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(D_elev),' deg )'],'fontsize',14);
axis([0 15 h_lims(1) h_lims(2)]);
grid on

subplot(322)
plot(time_axis,U_h(:,1));
hold on
plot(time_axis,U_h(:,2));
legend('Left Ear','Right Ear');
xlabel('Time [ms]','fontsize',12); ylabel('Amplitude','fontsize',12);
title(['IR amplitude ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(U_elev),' deg )'],'fontsize',14);
axis([0 15 h_lims(1) h_lims(2)]);
grid on

subplot(323)
plot(freq_axis,D_TF_mag_L);
hold on
plot(freq_axis,D_TF_mag_R);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Magnitude [dB]','fontsize',12);
title(['TF magnitude ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(D_elev),' deg )'],'fontsize',14);
axis([0 freq_axis(end) m_lims(1) m_lims(2)]);
grid on

subplot(324)
plot(freq_axis,U_TF_mag_L);
hold on
plot(freq_axis,U_TF_mag_R);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Magnitude [dB]','fontsize',12);
title(['TF magnitude ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(U_elev),' deg )'],'fontsize',14);
axis([0 freq_axis(end) m_lims(1) m_lims(2)]);
grid on

subplot(325)
plot(freq_axis,D_TF_pha_L);
hold on
plot(freq_axis,D_TF_pha_R);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Phase/\pi [rad/rad]','fontsize',12);
title(['TF phase ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(D_elev),' deg )'],'fontsize',14);
axis([0 freq_axis(end) p_lims(1) p_lims(2)]);
grid on

subplot(326)
plot(freq_axis,U_TF_pha_L);
hold on
plot(freq_axis,U_TF_pha_R);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Phase/\pi [rad/rad]','fontsize',12);
title(['TF phase ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(U_elev),' deg )'],'fontsize',14);
axis([0 freq_axis(end) p_lims(1) p_lims(2)]);
grid on
set(gcf,'position',[100 300 1200 1000]);

%% azimuth -180 ~ +180 deg, elevation -40 ~ +90 deg
% elevation -40 ~ +90 deg
%for elev = -40:5:90
for elev = U_elev
    fprintf('processing elevation %d deg\n',elev);
    
    IR_L = zeros(N_azim,N);
    IR_R = zeros(N_azim,N);
    TF_mag_L = zeros(N_azim,N/2);
    TF_mag_R = zeros(N_azim,N/2);
    TF_pha_L = zeros(N_azim,N/2);
    TF_pha_R = zeros(N_azim,N/2);
    
    % azimuth -180 ~ +180 deg
    for azim_idx = 1:N_azim
        azim = azims(azim_idx);
        
        % load IR
        [h,~] = ir_hats_F(azim,elev,root);
        IR_L(azim_idx,:) = h(:,1);
        IR_R(azim_idx,:) = h(:,2);
        
        % FFT
        H = fft(h);
        
        % magnitude
        TF_mag_L(azim_idx,:) = 20*log10(abs(H(1:N/2,1))/20e-6);
        TF_mag_R(azim_idx,:) = 20*log10(abs(H(1:N/2,2))/20e-6);
        
        % phase
        if pha_flag == 0
            TF_pha_L(azim_idx,:) = angle(H(1:N/2,1))/pi;
            TF_pha_R(azim_idx,:) = angle(H(1:N/2,2))/pi;
        else
            TF_pha_L(azim_idx,:) = unwrap(angle(H(1:N/2,1)))/pi;
            TF_pha_R(azim_idx,:) = unwrap(angle(H(1:N/2,2)))/pi;
            p_lims = [min([min(TF_pha_L), min(TF_pha_R)]), 0];
        end
    end
    
    % plot
    figure
    subplot(321)
    imagesc(time_axis,azims,IR_L,h_lims); axis xy
    colormap jet
    xlabel('Time [ms]','fontsize',12); ylabel('Azimuth [deg.]','fontsize',12);
    title(['Left ear IR amplitude ( Elevation ',num2str(elev),' deg )'],'fontsize',14);
    axis([0 15 azims(1) azims(end)]);
    grid on
    
    subplot(322)
    imagesc(time_axis,azims,IR_R,h_lims); axis xy
    colormap jet
    xlabel('Time [ms]','fontsize',12); ylabel('Azimuth [deg.]','fontsize',12);
    title(['Right ear IR amplitude ( Elevation ',num2str(elev),' deg )'],'fontsize',14);
    axis([0 15 azims(1) azims(end)]);
    grid on
    
    subplot(323)
    imagesc(freq_axis,azims,TF_mag_L,m_lims); axis xy
    colormap jet
    xlabel('Frequency [kHz]','fontsize',12); ylabel('Azimuth [deg.]','fontsize',12);
    title(['Left ear TF magnitude ( Elevation ',num2str(elev),' deg )'],'fontsize',14);
    axis([0 freq_axis(end) azims(1) azims(end)]);
    grid on
    
    subplot(324)
    imagesc(freq_axis,azims,TF_mag_R,m_lims); axis xy
    colormap jet
    xlabel('Frequency [kHz]','fontsize',12); ylabel('Azimuth [deg.]','fontsize',12);
    title(['Right ear TF magnitude ( Elevation ',num2str(elev),' deg )'],'fontsize',14);
    axis([0 freq_axis(end) azims(1) azims(end)]);
    grid on
    
    subplot(325)
    imagesc(freq_axis,azims,TF_pha_L,p_lims); axis xy
    colormap jet
    xlabel('Frequency [kHz]','fontsize',12); ylabel('Azimuth [deg.]','fontsize',12);
    title(['Left ear TF phase ( Elevation ',num2str(elev),' deg )'],'fontsize',14);
    axis([0 freq_axis(end) azims(1) azims(end)]);
    grid on
    
    subplot(326)
    imagesc(freq_axis,azims,TF_pha_R,p_lims); axis xy
    colormap jet
    xlabel('Frequency [kHz]','fontsize',12); ylabel('Azimuth [deg.]','fontsize',12);
    title(['Right ear TF phase ( Elevation ',num2str(elev),' deg )'],'fontsize',14);
    axis([0 freq_axis(end) azims(1) azims(end)]);
    grid on
    set(gcf,'position',[1300 300 1200 1000]);
end
