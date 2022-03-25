%-------------------------------------------------------------------------
%   Date : June 17, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : TF at reference point according to elevation
%   Synopsis : plot IR & TF at ref. point measured by Mic-R in 5 directions
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
h_lims = [-0.1 0.1];            % amplitude
m_lims = [60 100];              % magnitude [dB]
p_lims = [-1 1];                % phase/pi [rad/rad]

%% elevation X deg
% elevation
elev = 45;                      % -40 ~ +90 deg

% root of directory
root = '../00_Data/TF_Ref/Mic-R_5dir';

% load IR
[h,Fs] = ir_ref5_F(elev,root);
N = length(h);

% align at 3.25 ms (157th)
[h] = ir_align_F(h,157);

% time & freq. axis
time_axis = 1000*(0:N-1)/Fs;            % ms
freq_axis = 0.001*(0:N/2-1)*Fs/N;       % kHz

% magnitude
H = fft(h);
TF_mag = 20*log10(abs(H(1:N/2))/20e-6);

% phase
if pha_flag == 0                        % rad/rad
    TF_pha = angle(H(1:N/2))/pi;
else
    TF_pha = unwrap(angle(H(1:N/2)))/pi;
    p_lims = [min(TF_pha), 0];
end

% plot
figure
subplot(311)
plot(time_axis,h);
xlabel('Time [ms]','fontsize',12); ylabel('Amplitude','fontsize',12);
title(['IR amplitude at Ref. point ( Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([0 15 h_lims(1) h_lims(2)]);
grid on

subplot(312)
plot(freq_axis,TF_mag);
xlabel('Frequency [kHz]','fontsize',12); ylabel('Magnitude [dB]','fontsize',12);
title(['TF magnitude at Ref. point ( Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([0 freq_axis(end) m_lims(1) m_lims(2)]);
grid on

subplot(313)
plot(freq_axis,TF_pha);
xlabel('Frequency [kHz]','fontsize',12); ylabel('Phase/\pi [rad/rad]','fontsize',12);
title(['TF phase at Ref. point ( Elevation ',num2str(elev),' deg )'],'fontsize',14);
axis([0 freq_axis(end) p_lims(1) p_lims(2)]);
grid on
set(gcf,'position',[300 300 800 1000]);

%% elevation -40 ~ 90 deg
IR = zeros(N_elev,N);
TF_mag = zeros(N_elev,N/2);
TF_pha = zeros(N_elev,N/2);

for elev_idx = 1:N_elev
    % elevation
    elev = elevs(elev_idx);
    
    % load IR
    [h,~] = ir_ref5_F(elev,root);
    
    % align at 3.25 ms (157th)
    [h] = ir_align_F(h,157);
    IR(elev_idx,:) = h;
    
    % magnitude
    H = fft(h);
    TF_mag(elev_idx,:) = 20*log10(abs(H(1:N/2))/20e-6);
    
    % phase
    if pha_flag == 0
        TF_pha(elev_idx,:) = angle(H(1:N/2))/pi;
    else
        TF_pha(elev_idx,:) = unwrap(angle(H(1:N/2)))/pi;
        p_lims = [min(min(TF_pha)), 0];
    end
end

% plot
figure
subplot(311)
imagesc(time_axis,elevs,IR,h_lims); axis xy
colormap jet
xlabel('Time [ms]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
title('IR amplitude at Ref. point','fontsize',14);
axis([0 15 elevs(1) elevs(end)]);
grid on

subplot(312)
imagesc(freq_axis,elevs,TF_mag,m_lims); axis xy
colormap jet
xlabel('Frequency [kHz]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
title('TF magnitude at Ref. point','fontsize',14);
axis([0 freq_axis(end) elevs(1) elevs(end)]);
grid on

subplot(313)
imagesc(freq_axis,elevs,TF_pha,p_lims); axis xy
colormap jet
xlabel('Frequency [kHz]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
title('TF phase at Ref. point','fontsize',14);
axis([0 freq_axis(end) elevs(1) elevs(end)]);
grid on
set(gcf,'position',[1300 300 800 1000]);
