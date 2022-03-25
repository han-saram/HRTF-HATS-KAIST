%-------------------------------------------------------------------------
%   Date : July 06, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : HATS HRTF database generation
%   Synopsis : store 1 ms delayed HATS HRIR data in two columns of text file
%              such that left is first column, right is second column
%	Algorithm : -
%-------------------------------------------------------------------------

clc
clear
close all

% root of directory
root_hrtf = '../00_Data/HRTF_HATS';
root_hats = '../00_Data/TF_HATS';
root_ref5 = '../00_Data/TF_Ref/Mic-R_5dir';

% azimuth angles (0 ~ 355 deg)
for azim = 0:5:355
    fprintf('generating HATS HRTF database at %d deg azimuth...\n',azim);
    
    % make path
    path = sprintf('%s/a%03d', root_hrtf,azim);
    mkdir(path);
    
    % elevation angles (-40 ~ +90 deg)
    for elev = -40:5:90
        % load HATS IR
        [g_LR,~] = ir_hats_F(azim,elev,root_hats);
        g_L = g_LR(:,1);
        g_R = g_LR(:,2);
        
        % load Ref-point IR
        [g_0,~] = ir_ref5_F(elev,root_ref5);
        [g_0] = ir_align_F(g_0,157);        % align at 3.25 ms (157th)
        
        % windowing & zero-padding
        [g_L] = ir_window_F(g_L);
        [g_R] = ir_window_F(g_R);
        [g_0] = ir_window_F(g_0);
        
        % FFT
        G_L = fft(g_L);
        G_R = fft(g_R);
        G_0 = fft(g_0);
        
        % HRTF
        H_L = G_L./G_0;
        H_R = G_R./G_0;
        
        % HRIR
        h_L = ifft(H_L);
        h_R = ifft(H_R);
        
        % 1 ms delay for causality
        h_L = circshift(h_L,48);
        h_R = circshift(h_R,48);
        
        % pathname
        if elev >= 0
            name = sprintf('a%03de+%02d.txt', azim,elev);
        else
            name = sprintf('a%03de-%02d.txt', azim,abs(elev));
        end
        pathname = [path '/' name];
        
        % open file
        fid = fopen(pathname,'w');
        
        % write data
        for i = 1:length(h_R)
            fprintf(fid,'%+.18f %+.18f\n', h_L(i), h_R(i));
        end
        
        % close file
        fclose(fid);
    end
end
