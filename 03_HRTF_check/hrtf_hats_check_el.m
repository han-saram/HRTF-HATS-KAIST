%-------------------------------------------------------------------------
%   Date : July 06, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : HATS HRTF according to elevation
%   Synopsis : plot time & freq. responses of HATS HRTF measurements
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
h_lims = [-0.5 0.5];            % amplitude
m_lims = [-50 20];              % magnitude [dB]
p_lims = [-1 1];                % phase/pi [rad/rad]

%% azimuth 0 deg, elevation -X & X deg (median plane)
% source direction
azim = 0;
elev = 10;              % 0 ~ 40 deg

% load HRIR
U_elev = elev;
D_elev = -elev;
[U_h_L,U_h_R] = hrir_hats_F(azim,U_elev);
[D_h_L,D_h_R] = hrir_hats_F(azim,D_elev);
N = length(D_h_R);
Fs = 48e3;

% time & freq. axis
time_axis = 1000*(0:N-1)/Fs;            % ms
freq_axis = 0.001*(0:N/2-1)*Fs/N;       % kHz

% FFT
U_H_L = fft(U_h_L);
U_H_R = fft(U_h_R);
D_H_L = fft(D_h_L);
D_H_R = fft(D_h_R);

% magnitude
U_HRTF_mag_L = 20*log10(abs(U_H_L(1:N/2)));
U_HRTF_mag_R = 20*log10(abs(U_H_R(1:N/2)));
D_HRTF_mag_L = 20*log10(abs(D_H_L(1:N/2)));
D_HRTF_mag_R = 20*log10(abs(D_H_R(1:N/2)));

% phase
if pha_flag == 0                        % rad/rad
    U_HRTF_pha_L = angle(U_H_L(1:N/2))/pi;
    U_HRTF_pha_R = angle(U_H_R(1:N/2))/pi;
    D_HRTF_pha_L = angle(D_H_L(1:N/2))/pi;
    D_HRTF_pha_R = angle(D_H_R(1:N/2))/pi;
else
    U_HRTF_pha_L = unwrap(angle(U_H_L(1:N/2)))/pi;
    U_HRTF_pha_R = unwrap(angle(U_H_R(1:N/2)))/pi;
    D_HRTF_pha_L = unwrap(angle(D_H_L(1:N/2)))/pi;
    D_HRTF_pha_R = unwrap(angle(D_H_R(1:N/2)))/pi;
    p_lims = [min([min(U_HRTF_pha_L), min(U_HRTF_pha_R), min(D_HRTF_pha_L), min(D_HRTF_pha_R)]), 0];
end

% plot
figure
subplot(321)
plot(time_axis,D_h_L);
hold on
plot(time_axis,D_h_R);
legend('Left Ear','Right Ear');
xlabel('Time [ms]','fontsize',12); ylabel('Amplitude','fontsize',12);
title(['HRIR amplitude ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(D_elev),' deg )'],'fontsize',14);
axis([time_axis(1) time_axis(end) -1.0 1.0]);
grid on

subplot(322)
plot(time_axis,U_h_L);
hold on
plot(time_axis,U_h_R);
legend('Left Ear','Right Ear');
xlabel('Time [ms]','fontsize',12); ylabel('Amplitude','fontsize',12);
title(['HRIR amplitude ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(U_elev),' deg )'],'fontsize',14);
axis([time_axis(1) time_axis(end) -1.0 1.0]);
grid on

subplot(323)
plot(freq_axis,D_HRTF_mag_L);
hold on
plot(freq_axis,D_HRTF_mag_R);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Magnitude [dB]','fontsize',12);
title(['HRTF magnitude ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(D_elev),' deg )'],'fontsize',14);
axis([freq_axis(1) freq_axis(end) m_lims(1) m_lims(2)]);
grid on

subplot(324)
plot(freq_axis,U_HRTF_mag_L);
hold on
plot(freq_axis,U_HRTF_mag_R);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Magnitude [dB]','fontsize',12);
title(['HRTF magnitude ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(U_elev),' deg )'],'fontsize',14);
axis([freq_axis(1) freq_axis(end) m_lims(1) m_lims(2)]);
grid on

subplot(325)
plot(freq_axis,D_HRTF_pha_L);
hold on
plot(freq_axis,D_HRTF_pha_R);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Phase/\pi [rad/rad]','fontsize',12);
title(['HRTF phase ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(D_elev),' deg )'],'fontsize',14);
axis([freq_axis(1) freq_axis(end) p_lims(1) p_lims(2)]);
grid on

subplot(326)
plot(freq_axis,U_HRTF_pha_L);
hold on
plot(freq_axis,U_HRTF_pha_R);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Phase/\pi [rad/rad]','fontsize',12);
title(['HRTF phase ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(U_elev),' deg )'],'fontsize',14);
axis([freq_axis(1) freq_axis(end) p_lims(1) p_lims(2)]);
grid on
set(gcf,'position',[100 300 1200 1000]);

%% azimuth -180 ~ +180 deg, elevation -40 ~ +90 deg
% elevation -40 ~ +90 deg
for elev = -40:5:90
%for elev = U_elev
    fprintf('processing elevation %d deg\n',elev);
    
    HRIR_L = zeros(N_azim,N);
    HRIR_R = zeros(N_azim,N);
    HRTF_mag_L = zeros(N_azim,N/2);
    HRTF_mag_R = zeros(N_azim,N/2);
    HRTF_pha_L = zeros(N_azim,N/2);
    HRTF_pha_R = zeros(N_azim,N/2);
    
    % azimuth -180 ~ +180 deg
    for azim_idx = 1:N_azim
        azim = azims(azim_idx);
        
        % load HRIR
        [h_L,h_R] = hrir_hats_F(azim,elev);
        HRIR_L(azim_idx,:) = h_L;
        HRIR_R(azim_idx,:) = h_R;
        
        % FFT
        H_L = fft(h_L);
        H_R = fft(h_R);
        
        % magnitude
        HRTF_mag_L(azim_idx,:) = 20*log10(abs(H_L(1:N/2)));
        HRTF_mag_R(azim_idx,:) = 20*log10(abs(H_R(1:N/2)));
        
        % phase
        if pha_flag == 0
            HRTF_pha_L(azim_idx,:) = angle(H_L(1:N/2))/pi;
            HRTF_pha_R(azim_idx,:) = angle(H_R(1:N/2))/pi;
        else
            HRTF_pha_L(azim_idx,:) = unwrap(angle(H_L(1:N/2)))/pi;
            HRTF_pha_R(azim_idx,:) = unwrap(angle(H_R(1:N/2)))/pi;
            p_lims = [min([min(HRTF_pha_L), min(HRTF_pha_R)]), 0];
        end
    end
    
    % plot
    figure
    subplot(321)
    imagesc(time_axis,azims,HRIR_L,h_lims); axis xy
    colormap jet
    xlabel('Time [ms]','fontsize',12); ylabel('Azimuth [deg.]','fontsize',12);
    title(['Left ear HRIR amplitude ( Elevation ',num2str(elev),' deg )'],'fontsize',14);
    axis([0 4 azims(1) azims(end)]);
    grid on
    
    subplot(322)
    imagesc(time_axis,azims,HRIR_R,h_lims); axis xy
    colormap jet
    xlabel('Time [ms]','fontsize',12); ylabel('Azimuth [deg.]','fontsize',12);
    title(['Right ear HRIR amplitude ( Elevation ',num2str(elev),' deg )'],'fontsize',14);
    axis([0 4 azims(1) azims(end)]);
    grid on
    
    subplot(323)
    imagesc(freq_axis,azims,HRTF_mag_L,m_lims); axis xy
    colormap jet
    xlabel('Frequency [kHz]','fontsize',12); ylabel('Azimuth [deg.]','fontsize',12);
    title(['Left ear HRTF magnitude ( Elevation ',num2str(elev),' deg )'],'fontsize',14);
    axis([freq_axis(1) freq_axis(end) azims(1) azims(end)]);
    grid on
    
    subplot(324)
    imagesc(freq_axis,azims,HRTF_mag_R,m_lims); axis xy
    colormap jet
    xlabel('Frequency [kHz]','fontsize',12); ylabel('Azimuth [deg.]','fontsize',12);
    title(['Right ear HRTF magnitude ( Elevation ',num2str(elev),' deg )'],'fontsize',14);
    axis([freq_axis(1) freq_axis(end) azims(1) azims(end)]);
    grid on
    
    subplot(325)
    imagesc(freq_axis,azims,HRTF_pha_L,p_lims); axis xy
    colormap jet
    xlabel('Frequency [kHz]','fontsize',12); ylabel('Azimuth [deg.]','fontsize',12);
    title(['Left ear HRTF phase ( Elevation ',num2str(elev),' deg )'],'fontsize',14);
    axis([freq_axis(1) freq_axis(end) azims(1) azims(end)]);
    grid on
    
    subplot(326)
    imagesc(freq_axis,azims,HRTF_pha_R,p_lims); axis xy
    colormap jet
    xlabel('Frequency [kHz]','fontsize',12); ylabel('Azimuth [deg.]','fontsize',12);
    title(['Right ear HRTF phase ( Elevation ',num2str(elev),' deg )'],'fontsize',14);
    axis([freq_axis(1) freq_axis(end) azims(1) azims(end)]);
    grid on
    set(gcf,'position',[1300 300 1200 1000]);
end
