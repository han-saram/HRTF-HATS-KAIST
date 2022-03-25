%-------------------------------------------------------------------------
%   Date : July 16, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : HATS HRTF ILD calculation
%   Synopsis : calculate ILDs from HATS HRTF and store them in text files
%	Algorithm : -
%-------------------------------------------------------------------------

clc
clear
close all

% root path
root = '../00_Data/ILD_HATS';

% azimuth angles (0 ~ 355 deg)
for azim = 0:5:355
    fprintf('calculating ILDs at azimuth %d deg\n',azim);
    
    % make path
    path = sprintf('%s/a%03d', root,azim);
    mkdir(path);
    
    % elevation angles (-40 ~ +90 deg)
    for elev = -40:5:90
        % HRIR retrieval
        [h_L0,h_R0] = hrir_hats_F(azim,elev);
        N0 = length(h_R0);
        
        % zero padding
        N = 4800;
        h_L = zeros(N,1);
        h_R = zeros(N,1);
        h_L(1:N0) = h_L0;
        h_R(1:N0) = h_R0;
        
        % FFT
        H_L = fft(h_L);
        H_R = fft(h_R);
        
        % ILD
        ILD = 20*log10(abs(H_R(1:N/2)./H_L(1:N/2)));
        
        % freq. axis
        Fs = 48e3;
        f = (0:N/2-1)*Fs/N;
        
        % pathname
        if elev < 0
            name = sprintf('a%03de-%02d.txt', azim,abs(elev));
        else
            name = sprintf('a%03de+%02d.txt', azim,elev);
        end
        pathname = [path '/' name];
        
        % open file
        fid = fopen(pathname,'w');
        
        % write data
        for i = 1:N/2
            fprintf(fid,'%05d %+.3f\n', f(i), ILD(i));
        end
        
        % close file
        fclose(fid);
    end
end
