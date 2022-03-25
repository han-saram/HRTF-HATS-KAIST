%-------------------------------------------------------------------------
%   Date : July 06, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : HATS HRTF according to azimuth
%   Synopsis : plot time & freq. responses of HATS HRTF measurements
%	Algorithm : -
%-------------------------------------------------------------------------

clc
clear
close all

% elevation angles (-40 ~ +90 deg)
elevs = (-40:5:90)';
N_elev = length(elevs);

% phase flag (0: wrap, 1: unwrap)
pha_flag = 1;

% plot limit
h_lims = [-0.5 0.5];            % amplitude
m_lims = [-50 20];              % magnitude [dB]
p_lims = [-1 1];                % phase/pi [rad/rad]

%% azimuth -X & X deg, elevation 0 deg (horizontal plane)
% source direction
azim = 30;              % 0 ~ 180 deg
elev = 0;

% load HRIR
L_azim = -azim;
R_azim = azim;
[L_h_L,L_h_R] = hrir_hats_F(L_azim,elev);
[R_h_L,R_h_R] = hrir_hats_F(R_azim,elev);
N = length(R_h_R);
Fs = 48e3;

% time & freq. axis
time_axis = 1000*(0:N-1)/Fs;            % ms
freq_axis = 0.001*(0:N/2-1)*Fs/N;       % kHz

% FFT
L_H_L = fft(L_h_L);
L_H_R = fft(L_h_R);
R_H_L = fft(R_h_L);
R_H_R = fft(R_h_R);

% magnitude
L_HRTF_mag_L = 20*log10(abs(L_H_L(1:N/2)));
L_HRTF_mag_R = 20*log10(abs(L_H_R(1:N/2)));
R_HRTF_mag_L = 20*log10(abs(R_H_L(1:N/2)));
R_HRTF_mag_R = 20*log10(abs(R_H_R(1:N/2)));

% phase
if pha_flag == 0                        % rad/rad
    L_HRTF_pha_L = angle(L_H_L(1:N/2))/pi;
    L_HRTF_pha_R = angle(L_H_R(1:N/2))/pi;
    R_HRTF_pha_L = angle(R_H_L(1:N/2))/pi;
    R_HRTF_pha_R = angle(R_H_R(1:N/2))/pi;
else
    L_HRTF_pha_L = unwrap(angle(L_H_L(1:N/2)))/pi;
    L_HRTF_pha_R = unwrap(angle(L_H_R(1:N/2)))/pi;
    R_HRTF_pha_L = unwrap(angle(R_H_L(1:N/2)))/pi;
    R_HRTF_pha_R = unwrap(angle(R_H_R(1:N/2)))/pi;
    p_lims = [min([min(L_HRTF_pha_L), min(L_HRTF_pha_R), min(R_HRTF_pha_L), min(R_HRTF_pha_R)]), 0];
end

% plot
figure
subplot(321)
plot(time_axis,L_h_L);
hold on
plot(time_axis,L_h_R);
legend('Left Ear','Right Ear');
xlabel('Time [ms]','fontsize',12); ylabel('Amplitude','fontsize',12);
title(['HRIR amplitude ( Azimuth ',num2str(L_azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([time_axis(1) time_axis(end) -1.0 1.0]);
grid on

subplot(322)
plot(time_axis,R_h_L);
hold on
plot(time_axis,R_h_R);
legend('Left Ear','Right Ear');
xlabel('Time [ms]','fontsize',12); ylabel('Amplitude','fontsize',12);
title(['HRIR amplitude ( Azimuth ',num2str(R_azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([time_axis(1) time_axis(end) -1.0 1.0]);
grid on

subplot(323)
plot(freq_axis,L_HRTF_mag_L);
hold on
plot(freq_axis,L_HRTF_mag_R);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Magnitude [dB]','fontsize',12);
title(['HRTF magnitude ( Azimuth ',num2str(L_azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([freq_axis(1) freq_axis(end) m_lims(1) m_lims(2)]);
grid on

subplot(324)
plot(freq_axis,R_HRTF_mag_L);
hold on
plot(freq_axis,R_HRTF_mag_R);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Magnitude [dB]','fontsize',12);
title(['HRTF magnitude ( Azimuth ',num2str(R_azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([freq_axis(1) freq_axis(end) m_lims(1) m_lims(2)]);
grid on

subplot(325)
plot(freq_axis,L_HRTF_pha_L);
hold on
plot(freq_axis,L_HRTF_pha_R);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Phase/\pi [rad/rad]','fontsize',12);
title(['HRTF phase ( Azimuth ',num2str(L_azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([freq_axis(1) freq_axis(end) p_lims(1) p_lims(2)]);
grid on

subplot(326)
plot(freq_axis,R_HRTF_pha_L);
hold on
plot(freq_axis,R_HRTF_pha_R);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Phase/\pi [rad/rad]','fontsize',12);
title(['HRTF phase ( Azimuth ',num2str(R_azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([freq_axis(1) freq_axis(end) p_lims(1) p_lims(2)]);
grid on
set(gcf,'position',[100 300 1200 1000]);

%% azimuth 0 ~ 355 deg, elevation -40 ~ +90 deg
% azimuth 0 ~ 355 deg
%for azim = 0:5:355
for azim = R_azim
    fprintf('processing azimuth %d deg\n',azim);
    
    HRIR_L = zeros(N_elev,N);
    HRIR_R = zeros(N_elev,N);
    HRTF_mag_L = zeros(N_elev,N/2);
    HRTF_mag_R = zeros(N_elev,N/2);
    HRTF_pha_L = zeros(N_elev,N/2);
    HRTF_pha_R = zeros(N_elev,N/2);
    
    % elevation -40 ~ +90 deg
    for elev_idx = 1:N_elev
        elev = elevs(elev_idx);
        
        % load HRIR
        [h_L,h_R] = hrir_hats_F(azim,elev);
        HRIR_L(elev_idx,:) = h_L;
        HRIR_R(elev_idx,:) = h_R;
        
        % FFT
        H_L = fft(h_L);
        H_R = fft(h_R);
        
        % magnitude
        HRTF_mag_L(elev_idx,:) = 20*log10(abs(H_L(1:N/2)));
        HRTF_mag_R(elev_idx,:) = 20*log10(abs(H_R(1:N/2)));
        
        % phase
        if pha_flag == 0
            HRTF_pha_L(elev_idx,:) = angle(H_L(1:N/2))/pi;
            HRTF_pha_R(elev_idx,:) = angle(H_R(1:N/2))/pi;
        else
            HRTF_pha_L(elev_idx,:) = unwrap(angle(H_L(1:N/2)))/pi;
            HRTF_pha_R(elev_idx,:) = unwrap(angle(H_R(1:N/2)))/pi;
            p_lims = [min([min(HRTF_pha_L), min(HRTF_pha_R)]), 0];
        end
    end
    
    % plot
    figure
    subplot(321)
    imagesc(time_axis,elevs,HRIR_L,h_lims); axis xy
    colormap jet
    xlabel('Time [ms]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
    title(['Left ear HRIR amplitude ( Azimuth ',num2str(azim),' deg )'],'fontsize',14);
    axis([0 4 elevs(1) elevs(end)]);
    grid on
    
    subplot(322)
    imagesc(time_axis,elevs,HRIR_R,h_lims); axis xy
    colormap jet
    xlabel('Time [ms]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
    title(['Right ear HRIR amplitude ( Azimuth ',num2str(azim),' deg )'],'fontsize',14);
    axis([0 4 elevs(1) elevs(end)]);
    grid on
    
    subplot(323)
    imagesc(freq_axis,elevs,HRTF_mag_L,m_lims); axis xy
    colormap jet
    xlabel('Frequency [kHz]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
    title(['Left ear HRTF magnitude ( Azimuth ',num2str(azim),' deg )'],'fontsize',14);
    axis([freq_axis(1) freq_axis(end) elevs(1) elevs(end)]);
    grid on
    
    subplot(324)
    imagesc(freq_axis,elevs,HRTF_mag_R,m_lims); axis xy
    colormap jet
    xlabel('Frequency [kHz]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
    title(['Right ear HRTF magnitude ( Azimuth ',num2str(azim),' deg )'],'fontsize',14);
    axis([freq_axis(1) freq_axis(end) elevs(1) elevs(end)]);
    grid on
    
    subplot(325)
    imagesc(freq_axis,elevs,HRTF_pha_L,p_lims); axis xy
    colormap jet
    xlabel('Frequency [kHz]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
    title(['Left ear HRTF phase ( Azimuth ',num2str(azim),' deg )'],'fontsize',14);
    axis([freq_axis(1) freq_axis(end) elevs(1) elevs(end)]);
    grid on
    
    subplot(326)
    imagesc(freq_axis,elevs,HRTF_pha_R,p_lims); axis xy
    colormap jet
    xlabel('Frequency [kHz]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
    title(['Right ear HRTF phase ( Azimuth ',num2str(azim),' deg )'],'fontsize',14);
    axis([freq_axis(1) freq_axis(end) elevs(1) elevs(end)]);
    grid on
    set(gcf,'position',[1300 300 1200 1000]);
end
