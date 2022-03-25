%-------------------------------------------------------------------------
%   Date : July 06, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------

clc
clear
close all

% phase flag (0: wrap, 1: unwrap)
pha_flag = 1;

% plot limit
h_lims = [-1.0 1.5];            % amplitude
m_lims = [-30 20];              % magnitude [dB]
p_lims = [-1 1];                % phase/pi [rad/rad]

%% raw HRIR vs 1 ms delayed HRIR
% source direction
azim = 90;
elev = 0;

% load raw HRIR
[h_L_raw,h_R_raw] = hrir_hats_F(azim,elev,0);
[h_L_1ms,h_R_1ms] = hrir_hats_F(azim,elev);
N = length(h_R_1ms);
Fs = 48e3;

% time & freq. axis
time_axis = 1000*(0:N-1)/Fs;            % ms
freq_axis = 0.001*(0:N/2-1)*Fs/N;       % kHz

% FFT
H_L_raw = fft(h_L_raw);
H_R_raw = fft(h_R_raw);
H_L_1ms = fft(h_L_1ms);
H_R_1ms = fft(h_R_1ms);

% magnitude
HRTF_mag_L_raw = 20*log10(abs(H_L_raw(1:N/2)));
HRTF_mag_R_raw = 20*log10(abs(H_R_raw(1:N/2)));
HRTF_mag_L_1ms = 20*log10(abs(H_L_1ms(1:N/2)));
HRTF_mag_R_1ms = 20*log10(abs(H_R_1ms(1:N/2)));

% phase
if pha_flag == 0                        % rad/rad
    HRTF_pha_L_raw = angle(H_L_raw(1:N/2))/pi;
    HRTF_pha_R_raw = angle(H_R_raw(1:N/2))/pi;
    HRTF_pha_L_1ms = angle(H_L_1ms(1:N/2))/pi;
    HRTF_pha_R_1ms = angle(H_R_1ms(1:N/2))/pi;
else
    HRTF_pha_L_raw = unwrap(angle(H_L_raw(1:N/2)))/pi;
    HRTF_pha_R_raw = unwrap(angle(H_R_raw(1:N/2)))/pi;
    HRTF_pha_L_1ms = unwrap(angle(H_L_1ms(1:N/2)))/pi;
    HRTF_pha_R_1ms = unwrap(angle(H_R_1ms(1:N/2)))/pi;
    p_lims = [-40, 10];
end

% plot
figure
subplot(311)
plot(time_axis,h_L_raw);
hold on
plot(time_axis,h_R_raw);
legend('Left Ear','Right Ear');
xlabel('Time [ms]','fontsize',12); ylabel('Amplitude','fontsize',12);
title(['Raw HRIR amplitude ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([time_axis(1) time_axis(end) h_lims(1) h_lims(2)]);
grid on

subplot(312)
semilogx(freq_axis,HRTF_mag_L_raw);
hold on
semilogx(freq_axis,HRTF_mag_R_raw);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Magnitude [dB]','fontsize',12);
title(['Raw HRTF magnitude ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([0 20 m_lims(1) m_lims(2)]);
grid on

subplot(313)
semilogx(freq_axis,HRTF_pha_L_raw);
hold on
semilogx(freq_axis,HRTF_pha_R_raw);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Phase/\pi [rad/rad]','fontsize',12);
title(['Raw HRTF phase ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([0 20 p_lims(1) p_lims(2)]);
grid on
set(gcf,'position',[100 300 700 1000]);

figure
subplot(311)
plot(time_axis,h_L_1ms);
hold on
plot(time_axis,h_R_1ms);
legend('Left Ear','Right Ear');
xlabel('Time [ms]','fontsize',12); ylabel('Amplitude','fontsize',12);
title(['Delayed HRIR amplitude ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([time_axis(1) time_axis(end) h_lims(1) h_lims(2)]);
grid on

subplot(312)
semilogx(freq_axis,HRTF_mag_L_1ms);
hold on
semilogx(freq_axis,HRTF_mag_R_1ms);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Magnitude [dB]','fontsize',12);
title(['Delayed HRTF magnitude ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([0 20 m_lims(1) m_lims(2)]);
grid on

subplot(313)
semilogx(freq_axis,HRTF_pha_L_1ms);
hold on
semilogx(freq_axis,HRTF_pha_R_1ms);
legend('Left Ear','Right Ear','Location','southwest');
xlabel('Frequency [kHz]','fontsize',12); ylabel('Phase/\pi [rad/rad]','fontsize',12);
title(['Delayed HRTF phase ( Azimuth ',num2str(azim),' deg, Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([0 20 p_lims(1) p_lims(2)]);
grid on
set(gcf,'position',[800 300 700 1000]);

%% raw HRIR vs 1 ms delayed HRIR - convolution of square pulse
% square pulse
x = zeros(100*48,1);            % length : 100 ms
x(30*48+1:60*48) = 1;           % width : 30 ~ 60 ms
N = length(x);

% time & freq. axis
t = 1000*(0:N-1)/Fs;            % ms
f = 0.001*(0:N/2-1)*Fs/N;       % kHz

% convolution
if (azim > 180) || (azim < 0)
    x_raw = conv(x,h_L_raw);
    x_1ms = conv(x,h_L_1ms);
else
    x_raw = conv(x,h_R_raw);
    x_1ms = conv(x,h_R_1ms);
end
x_raw = x_raw(1:N);
x_1ms = x_1ms(1:N);

% FFT
X = fft(x);
X_raw = fft(x_raw);
X_1ms = fft(x_1ms);

% magnitude
HRTF_mag = 20*log10(abs(X(1:N/2)));
HRTF_mag_raw = 20*log10(abs(X_raw(1:N/2)));
HRTF_mag_1ms = 20*log10(abs(X_1ms(1:N/2)));

% phase
if pha_flag == 0                        % rad/rad
    HRTF_pha = angle(X(1:N/2))/pi;
    HRTF_pha_raw = angle(X_raw(1:N/2))/pi;
    HRTF_pha_1ms = angle(X_1ms(1:N/2))/pi;
else
    HRTF_pha = unwrap(angle(X(1:N/2)))/pi;
    HRTF_pha_raw = unwrap(angle(X_raw(1:N/2)))/pi;
    HRTF_pha_1ms = unwrap(angle(X_1ms(1:N/2)))/pi;
    p_lim1 = min([min(HRTF_pha), min(HRTF_pha_raw), min(HRTF_pha_1ms)]);
    p_lim2 = max([max(HRTF_pha), max(HRTF_pha_raw), max(HRTF_pha_1ms)]);
    p_lims = [p_lim1, p_lim2];
end

% plot
figure
subplot(331)
plot(t,x);
xlabel('Time [ms]','fontsize',12); ylabel('Amplitude','fontsize',12);
title('Square pulse','fontsize',14);
axis([t(1) t(end) -3 3]);
grid on

subplot(332)
plot(t,x_raw);
xlabel('Time [ms]','fontsize',12); ylabel('Amplitude','fontsize',12);
title('Square pulse convolved with raw ipsilateral HRIR','fontsize',14);
axis([t(1) t(end) -3 3]);
grid on

subplot(333)
plot(t,x_1ms);
xlabel('Time [ms]','fontsize',12); ylabel('Amplitude','fontsize',12);
title('Square pulse convolved with delayed ipsilateral HRIR','fontsize',14);
axis([t(1) t(end) -3 3]);
grid on

subplot(334)
semilogx(f,HRTF_mag);
xlabel('Frequency [kHz]','fontsize',12); ylabel('Magnitude [dB]','fontsize',12);
axis([0 20 -20 60]);
grid on

subplot(335)
semilogx(f,HRTF_mag_raw);
xlabel('Frequency [kHz]','fontsize',12); ylabel('Magnitude [dB]','fontsize',12);
axis([0 20 -20 60]);
grid on

subplot(336)
semilogx(f,HRTF_mag_1ms);
xlabel('Frequency [kHz]','fontsize',12); ylabel('Magnitude [dB]','fontsize',12);
axis([0 20 -20 60]);
grid on

subplot(337)
semilogx(f,HRTF_pha);
xlabel('Frequency [kHz]','fontsize',12); ylabel('Phase/\pi [rad/rad]','fontsize',12);
axis([0 20 p_lims(1) p_lims(2)]);
grid on

subplot(338)
semilogx(f,HRTF_pha_raw);
xlabel('Frequency [kHz]','fontsize',12); ylabel('Phase/\pi [rad/rad]','fontsize',12);
axis([0 20 p_lims(1) p_lims(2)]);
grid on

subplot(339)
semilogx(f,HRTF_pha_1ms);
xlabel('Frequency [kHz]','fontsize',12); ylabel('Phase/\pi [rad/rad]','fontsize',12);
axis([0 20 p_lims(1) p_lims(2)]);
grid on
set(gcf,'position',[300 300 2000 1000]);
