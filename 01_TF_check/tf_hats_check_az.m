%-------------------------------------------------------------------------
%   Date : June 03, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : HATS TF according to azimuth
%   Synopsis : plot time & freq. responses of HATS TF measurements
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
h_lims = [-0.05 0.05];          % amplitude
m_lims = [40 100];              % magnitude [dB]
p_lims = [-1 1];                % phase/pi [rad/rad]

%% azimuth -X & X deg, elevation 0 deg (horizontal plane)
% source direction
azim = 30;              % 0 ~ 180 deg
elev = 0;

% root of directory
root = '../00_Data/TF_HATS';

% load IR
L_azim = -azim;
R_azim = azim;
[L_h,Fs] = ir_hats_F(L_azim,elev,root);
[R_h,~] = ir_hats_F(R_azim,elev,root);

% length of IR
N = length(R_h);

% time & freq. axis
time_axis = 1000*(0:N-1)/Fs;            % ms
freq_axis = 0.001*(0:N/2-1)*Fs/N;       % kHz

% FFT
L_H = fft(L_h);
R_H = fft(R_h);

% magnitude
L_TF_mag_L = 20*log10(abs(L_H(1:N/2,1))/20e-6);
L_TF_mag_R = 20*log10(abs(L_H(1:N/2,2))/20e-6);
R_TF_mag_L = 20*log10(abs(R_H(1:N/2,1))/20e-6);
R_TF_mag_R = 20*log10(abs(R_H(1:N/2,2))/20e-6);

% phase
if pha_flag == 0                        % rad/rad
    L_TF_pha_L = angle(L_H(1:N/2,1))/pi;
    L_TF_pha_R = angle(L_H(1:N/2,2))/pi;
    R_TF_pha_L = angle(R_H(1:N/2,1))/pi;
    R_TF_pha_R = angle(R_H(1:N/2,2))/pi;
else
    L_TF_pha_L = unwrap(angle(L_H(1:N/2,1)))/pi;
    L_TF_pha_R = unwrap(angle(L_H(1:N/2,2)))/pi;
    R_TF_pha_L = unwrap(angle(R_H(1:N/2,1)))/pi;
    R_TF_pha_R = unwrap(angle(R_H(1:N/2,2)))/pi;
    p_lims = [min([min(L_TF_pha_L), min(L_TF_pha_R), min(R_TF_pha_L), min(R_TF_pha_R)]), 0];
end

% plot
figure
subplot(321)
plot(time_axis,L_h(:,1));
hold on
plot(time_axis,L_h(:,2));
legend('Left Ear','Right Ear');
xlabel('Time [ms]','fontsize',12); ylabel('Amplitude','fontsize',12);
title(['IR amplitude ( Azimuth ',num2str(L_azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([0 15 h_lims(1) h_lims(2)]);
grid on

subplot(322)
plot(time_axis,R_h(:,1));
hold on
plot(time_axis,R_h(:,2));
legend('Left Ear','Right Ear');
xlabel('Time [ms]','fontsize',12); ylabel('Amplitude','fontsize',12);
title(['IR amplitude ( Azimuth ',num2str(R_azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([0 15 h_lims(1) h_lims(2)]);
grid on

subplot(323)
plot(freq_axis,L_TF_mag_L);
hold on
plot(freq_axis,L_TF_mag_R);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Magnitude [dB]','fontsize',12);
title(['TF magnitude ( Azimuth ',num2str(L_azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([0 freq_axis(end) m_lims(1) m_lims(2)]);
grid on

subplot(324)
plot(freq_axis,R_TF_mag_L);
hold on
plot(freq_axis,R_TF_mag_R);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Magnitude [dB]','fontsize',12);
title(['TF magnitude ( Azimuth ',num2str(R_azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([0 freq_axis(end) m_lims(1) m_lims(2)]);
grid on

subplot(325)
plot(freq_axis,L_TF_pha_L);
hold on
plot(freq_axis,L_TF_pha_R);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Phase/\pi [rad/rad]','fontsize',12);
title(['TF phase ( Azimuth ',num2str(L_azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([0 freq_axis(end) p_lims(1) p_lims(2)]);
grid on

subplot(326)
plot(freq_axis,R_TF_pha_L);
hold on
plot(freq_axis,R_TF_pha_R);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Phase/\pi [rad/rad]','fontsize',12);
title(['TF phase ( Azimuth ',num2str(R_azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([0 freq_axis(end) p_lims(1) p_lims(2)]);
grid on
set(gcf,'position',[100 300 1200 1000]);

%% azimuth 0 ~ 355 deg, elevation -40 ~ +90 deg
% azimuth 0 ~ 355 deg
%for azim = 0:5:355
for azim = R_azim
    fprintf('processing azimuth %d deg\n',azim);
    
    IR_L = zeros(N_elev,N);
    IR_R = zeros(N_elev,N);
    TF_mag_L = zeros(N_elev,N/2);
    TF_mag_R = zeros(N_elev,N/2);
    TF_pha_L = zeros(N_elev,N/2);
    TF_pha_R = zeros(N_elev,N/2);
    
    % elevation -40 ~ +90 deg
    for elev_idx = 1:N_elev
        elev = elevs(elev_idx);
        
        % load IR
        [h,~] = ir_hats_F(azim,elev,root);
        IR_L(elev_idx,:) = h(:,1);
        IR_R(elev_idx,:) = h(:,2);
        
        % FFT
        H = fft(h);
        
        % magnitude
        TF_mag_L(elev_idx,:) = 20*log10(abs(H(1:N/2,1))/20e-6);
        TF_mag_R(elev_idx,:) = 20*log10(abs(H(1:N/2,2))/20e-6);
        
        % phase
        if pha_flag == 0
            TF_pha_L(elev_idx,:) = angle(H(1:N/2,1))/pi;
            TF_pha_R(elev_idx,:) = angle(H(1:N/2,2))/pi;
        else
            TF_pha_L(elev_idx,:) = unwrap(angle(H(1:N/2,1)))/pi;
            TF_pha_R(elev_idx,:) = unwrap(angle(H(1:N/2,2)))/pi;
            p_lims = [min([min(TF_pha_L), min(TF_pha_R)]), 0];
        end
    end
    
    % plot
    figure
    subplot(321)
    imagesc(time_axis,elevs,IR_L,h_lims); axis xy
    colormap jet
    xlabel('Time [ms]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
    title(['Left ear IR amplitude ( Azimuth ',num2str(azim),' deg )'],'fontsize',14);
    axis([0 15 elevs(1) elevs(end)]);
    grid on
    
    subplot(322)
    imagesc(time_axis,elevs,IR_R,h_lims); axis xy
    colormap jet
    xlabel('Time [ms]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
    title(['Right ear IR amplitude ( Azimuth ',num2str(azim),' deg )'],'fontsize',14);
    axis([0 15 elevs(1) elevs(end)]);
    grid on
    
    subplot(323)
    imagesc(freq_axis,elevs,TF_mag_L,m_lims); axis xy
    colormap jet
    xlabel('Frequency [kHz]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
    title(['Left ear TF magnitude ( Azimuth ',num2str(azim),' deg )'],'fontsize',14);
    axis([0 freq_axis(end) elevs(1) elevs(end)]);
    grid on
    
    subplot(324)
    imagesc(freq_axis,elevs,TF_mag_R,m_lims); axis xy
    colormap jet
    xlabel('Frequency [kHz]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
    title(['Right ear TF magnitude ( Azimuth ',num2str(azim),' deg )'],'fontsize',14);
    axis([0 freq_axis(end) elevs(1) elevs(end)]);
    grid on
    
    subplot(325)
    imagesc(freq_axis,elevs,TF_pha_L,p_lims); axis xy
    colormap jet
    xlabel('Frequency [kHz]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
    title(['Left ear TF phase ( Azimuth ',num2str(azim),' deg )'],'fontsize',14);
    axis([0 freq_axis(end) elevs(1) elevs(end)]);
    grid on
    
    subplot(326)
    imagesc(freq_axis,elevs,TF_pha_R,p_lims); axis xy
    colormap jet
    xlabel('Frequency [kHz]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
    title(['Right ear TF phase ( Azimuth ',num2str(azim),' deg )'],'fontsize',14);
    axis([0 freq_axis(end) elevs(1) elevs(end)]);
    grid on
    set(gcf,'position',[1300 300 1200 1000]);
end
