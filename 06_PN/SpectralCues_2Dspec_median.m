%-------------------------------------------------------------------------
%   Date : July 22, 2021
%   Place : Human Lab, KAIST
%   Coder : Ko, Byeong Yun & Lee, Gyeong Tae 
%-------------------------------------------------------------------------
%	Title : Spectral cues of HATS HRTF
%   Synopsis : plot spectral cues depending on azi & ele
%	Algorithm : -
%-------------------------------------------------------------------------

clc
clear
close all

% elevation angles (-40 ~ +90 deg)
elevs = (-40:10:90)';
N_elev = length(elevs);

% phase flag (0: wrap, 1: unwrap)
pha_flag = 1;

% plot limit
m_lims = [-40 20];              % magnitude [dB]

%% individual HRTF
% source direction
azim = 0;              % 0 ~ 180 deg
elev = 0;

% load HRIR
L_azim = -azim;
[L_h_L,L_h_R] = hrir_hats_F(L_azim,elev);
L_h_L = L_h_L./max(abs(L_h_L));
N = length(L_h_L);
Fs = 48e3;

origin_L_h_L = L_h_L;

% windowing
n = 96*1;
win = hann(n);
win = [win; zeros(N-n,1)];

L_h_L = L_h_L .* circshift(win,find(L_h_L == max(L_h_L))-n/2); 

% time & freq. axis
time_axis = 1000*(0:N-1)/Fs;            % ms
freq_axis = 0.001*(0:N/2-1)*Fs/N;       % kHz

% FFT
L_H_L = fft(L_h_L);
origin_L_H_L = fft(origin_L_h_L);

% magnitude
L_HRTF_mag_L = 20*log10(abs(L_H_L(1:N/2)));
origin_L_HRTF_mag_L = 20*log10(abs(origin_L_H_L(1:N/2)));

% spectral cues
P_L_HRTF_mag_L = islocalmax(L_HRTF_mag_L);
N_L_HRTF_mag_L = islocalmin(L_HRTF_mag_L);

P_origin_L_HRTF_mag_L = islocalmax(origin_L_HRTF_mag_L);
N_origin_L_HRTF_mag_L = islocalmin(origin_L_HRTF_mag_L);


temp = find(P_L_HRTF_mag_L == 1);
P_L_HRTF_mag_L(temp(find(freq_axis(P_L_HRTF_mag_L) < 3))) = 0;
temp = find(N_L_HRTF_mag_L == 1);
N_L_HRTF_mag_L(temp(find(freq_axis(N_L_HRTF_mag_L) < 3))) = 0;

temp = find(P_origin_L_HRTF_mag_L == 1);
P_origin_L_HRTF_mag_L(temp(find(freq_axis(P_origin_L_HRTF_mag_L) < 3))) = 0;
temp = find(N_origin_L_HRTF_mag_L == 1);
N_origin_L_HRTF_mag_L(temp(find(freq_axis(N_origin_L_HRTF_mag_L) < 3))) = 0;

% plot
figure
subplot(121)
plot(time_axis,origin_L_h_L,'k','LineWidth',1.8);
hold on
plot(time_axis,L_h_L,'g-.','LineWidth',1.8);
axis([time_axis(1) time_axis(end) -1.0 1.0]);
set(gca,'FontName','times','FontSize',18);
legend('Original','Windowed','FontSize',14);
xlabel('Time, ms','FontSize',18)
ylabel('Amplitude','FontSize',18)
title(['HRIR amplitude at \theta=',num2str(L_azim),'\circ, \phi=',num2str(elev),'\circ'],'FontSize',18)
grid on

subplot(122)
plot(freq_axis,origin_L_HRTF_mag_L,'k','LineWidth',1.8);
hold on
plot(freq_axis,L_HRTF_mag_L,'g-.','LineWidth',1.8);
hold on
plot(freq_axis(P_L_HRTF_mag_L),L_HRTF_mag_L(P_L_HRTF_mag_L),'x','color','#000000','MarkerSize', 13,'LineWidth', 1.8);
hold on
plot(freq_axis(N_L_HRTF_mag_L),L_HRTF_mag_L(N_L_HRTF_mag_L),'o','color','#000000','MarkerSize', 10,'LineWidth', 1.8);

axis([freq_axis(1) freq_axis(end) m_lims(1) m_lims(2)]);
grid on

set(gca,'FontName','times','FontSize',18);
legend('Original','Windowed','Peaks','Notches','Location','southwest','FontSize',14);
xlabel('Frequency, kHz','FontSize',18)
ylabel('Magnitude, dB','FontSize',18)
title(['HRTF spectrum at \theta=',num2str(L_azim),'\circ, \phi=',num2str(elev),'\circ'],'FontSize',18)
set(gcf,'position',[50 100 1200 350]);



%% median plane
R_azim = 0;

% windowing
n = 96*1;
win = hann(n);
win = [win; zeros(N-n,1)];

spec_cues = zeros(2*N_elev-1,12);
HRTF_L = zeros(2*N_elev-1,N/2);
HRTF_R = zeros(2*N_elev-1,N/2);

% median plane
azims = [0 180];
for i = 1:2
    azim = azims(i); 
    fprintf('processing median plane...');   
        
    % elevation -40 ~ +90 deg
    for elev_idx = 1:N_elev
        elev = elevs(elev_idx);
        
        % load HRIR
        [h_L,h_R] = hrir_hats_F(azim,elev);
        N = length(h_L);
        h_L = h_L./max(abs(h_L));
        h_R = h_R./max(abs(h_R));
        Fs = 48e3;

        h_L = h_L .* circshift(win,find(h_L == max(h_L))-n/2); 
        h_R = h_R .* circshift(win,find(h_R == max(h_R))-n/2); 
        
        % time & freq. axis
        time_axis = 1000*(0:N-1)/Fs;            % ms
        freq_axis = 0.001*(0:N/2-1)*Fs/N;       % kHz
                
        % FFT
        H_L = fft(h_L);
        H_R = fft(h_R);
        
        if i == 1
            % magnitude
            HRTF_L(elev_idx,:) = 20*log10(abs(H_L(1:N/2)));
            HRTF_R(elev_idx,:) = 20*log10(abs(H_R(1:N/2)));
        else
            % magnitude
            HRTF_L(2*N_elev-elev_idx,:) = 20*log10(abs(H_L(1:N/2)));
            HRTF_R(2*N_elev-elev_idx,:) = 20*log10(abs(H_R(1:N/2)));
        end
             

        % spectral cues
        if i == 1
            P_HRTF_L = islocalmax(HRTF_L(elev_idx,:));
            N_HRTF_L = islocalmin(HRTF_L(elev_idx,:));
            P_HRTF_R = islocalmax(HRTF_R(elev_idx,:));
            N_HRTF_R = islocalmin(HRTF_R(elev_idx,:));
        else
            P_HRTF_L = islocalmax(HRTF_L(2*N_elev-elev_idx,:));
            N_HRTF_L = islocalmin(HRTF_L(2*N_elev-elev_idx,:));
            P_HRTF_R = islocalmax(HRTF_R(2*N_elev-elev_idx,:));
            N_HRTF_R = islocalmin(HRTF_R(2*N_elev-elev_idx,:));
        end

        temp = find(P_HRTF_L == 1);
        P_HRTF_L(temp(find(freq_axis(P_HRTF_L) < 3 | freq_axis(P_HRTF_L) > 16.5))) = 0;
        temp = find(N_HRTF_L == 1);
        N_HRTF_L(temp(find(freq_axis(N_HRTF_L) < 3 | freq_axis(N_HRTF_L) > 16.5))) = 0;
        temp = find(P_HRTF_R == 1);
        P_HRTF_R(temp(find(freq_axis(P_HRTF_R) < 3 | freq_axis(P_HRTF_R) > 16.5))) = 0;
        temp = find(N_HRTF_R == 1);
        N_HRTF_R(temp(find(freq_axis(N_HRTF_R) < 3 | freq_axis(N_HRTF_R) > 16.5))) = 0;
        
        if i == 1
            temp = freq_axis(P_HRTF_L);
            if length(temp) == 2
                spec_cues(elev_idx,1:2) = temp(1:2);
                spec_cues(elev_idx,3) = temp(2);
            elseif length(temp) == 1
                spec_cues(elev_idx,1:3) = temp(1);
            else
                spec_cues(elev_idx,1:3) = temp(1:3);            
            end
            temp = freq_axis(N_HRTF_L);
            if length(temp) == 2
                spec_cues(elev_idx,4:5) = temp(1:2);
                spec_cues(elev_idx,6) = temp(2);
            elseif length(temp) == 1
                spec_cues(elev_idx,4:6) = temp(1);
            else
                spec_cues(elev_idx,4:6) = temp(1:3);
            end
            temp = freq_axis(P_HRTF_R);
            if length(temp) == 2
                spec_cues(elev_idx,7:8) = temp(1:2);
                spec_cues(elev_idx,9) = temp(2);
            elseif length(temp) == 1
                spec_cues(elev_idx,7:9) = temp(1);
            else
                spec_cues(elev_idx,7:9) = temp(1:3);
            end
            temp = freq_axis(N_HRTF_R);
            if length(temp) == 2
                spec_cues(elev_idx,10:11) = temp(1:2);
                spec_cues(elev_idx,12) = temp(2);
            elseif length(temp) == 1
                spec_cues(elev_idx,10:12) = temp(1);
            else
                spec_cues(elev_idx,10:12) = temp(1:3);
            end
            
        else
            temp = freq_axis(P_HRTF_L);
            if length(temp) == 2
                spec_cues(2*N_elev-elev_idx,1:2) = temp(1:2);
                spec_cues(2*N_elev-elev_idx,3) = temp(2);
            elseif length(temp) == 1
                spec_cues(2*N_elev-elev_idx,1:3) = temp(1);
            else
                spec_cues(2*N_elev-elev_idx,1:3) = temp(1:3);            
            end
            temp = freq_axis(N_HRTF_L);
            if length(temp) == 2
                spec_cues(2*N_elev-elev_idx,4:5) = temp(1:2);
                spec_cues(2*N_elev-elev_idx,6) = temp(2);
            elseif length(temp) == 1
                spec_cues(2*N_elev-elev_idx,4:6) = temp(1);
            else
                spec_cues(2*N_elev-elev_idx,4:6) = temp(1:3);
            end
            temp = freq_axis(P_HRTF_R);
            if length(temp) == 2
                spec_cues(2*N_elev-elev_idx,7:8) = temp(1:2);
                spec_cues(2*N_elev-elev_idx,9) = temp(2);
            elseif length(temp) == 1
                spec_cues(2*N_elev-elev_idx,7:9) = temp(1);
            else
                spec_cues(2*N_elev-elev_idx,7:9) = temp(1:3);
            end
            temp = freq_axis(N_HRTF_R);
            if length(temp) == 2
                spec_cues(2*N_elev-elev_idx,10:11) = temp(1:2);
                spec_cues(2*N_elev-elev_idx,12) = temp(2);
            elseif length(temp) == 1
                spec_cues(2*N_elev-elev_idx,10:12) = temp(1);
            elseif isempty(temp)
                
            else
                spec_cues(2*N_elev-elev_idx,10:12) = temp(1:3);
            end
        end
    end
end    

P1_L = spec_cues(:,1);
P2_L = spec_cues(:,2);
P3_L = spec_cues(:,3);
N1_L = spec_cues(:,4);
N2_L = spec_cues(:,5);
N3_L = spec_cues(:,6);
P1_R = spec_cues(:,7);
P2_R = spec_cues(:,8);
P3_R = spec_cues(:,9);
N1_R = spec_cues(:,10);
N2_R = spec_cues(:,11);
N3_R = spec_cues(:,12);
elevs = [-40 : 10 : 220]';

Ele = elevs;
T = table(Ele,P1_L,P2_L,P3_L,N1_L,N2_L,N3_L,P1_R,P2_R,P3_R,N1_R,N2_R,N3_R);
eval(['writetable(T,''..\00_Data\PN_HATS\HRTF_cues_median\median.txt'',''Delimiter'',''\t'',''WriteRowNames'',true);'])

Ang = reshape([elevs elevs elevs],[],1);
Peak = reshape(spec_cues(:,1:3),[],1);
Notch = reshape(spec_cues(:,4:6),[],1);
markersize1 = 100*ones(length(Peak),1);
markersize2 = 60*ones(length(Peak),1);
figure
im = imagesc(freq_axis,elevs,HRTF_L,m_lims); axis xy
im.AlphaData = .8;
colormap jet
hold on
scatter(Peak,Ang,markersize1,'x','MarkerEdgeColor','#000000','LineWidth', 1.8)
hold on
scatter(Notch,Ang,markersize2,'o','MarkerEdgeColor','#000000','LineWidth', 1.8)    

set(gca,'FontName','times','FontSize',18);
title(['Spectral cues of left ear at median plane'],'FontSize',18)
xlabel('Frequency, kHz','FontSize',18)
ylabel('\phi, deg','FontSize',18)
axis([freq_axis(1) freq_axis(end) elevs(1) elevs(end)]);
xlim([0 20])
set(gcf,'position',[50 100 650 400]);
colorbar
hc=colorbar;
title(hc,'dB');
grid on

