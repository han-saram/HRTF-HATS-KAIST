%-------------------------------------------------------------------------
%   Date : July 16, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : HATS HRTF full band ILD calculation
%   Synopsis : calculate ILDs from HATS HRTF and store them in text files
%	Algorithm : -
%-------------------------------------------------------------------------

clc
clear
close all

% make path
path = '../00_Data/ILD_HATS_full';
mkdir(path);

% azimuth angles (0 ~ 355 deg)
for azim = 0:5:355
    fprintf('calculating ILDs at azimuth %d deg\n',azim);
    
    % pathname
    name = sprintf('a%03d.txt', azim);
    pathname = [path '/' name];
    
    % open file
    fid = fopen(pathname,'w');
    
    % elevation angles (-40 ~ +90 deg)
    for elev = -40:5:90
        % HRIR retrieval
        [h_L,h_R] = hrir_hats_F(azim,elev);
        N = length(h_R);
        
        % FFT
        H_L = fft(h_L);
        H_R = fft(h_R);
        
        % power sum
        num = 0;
        den = 0;
        for k = 1:N/2
            num = num + abs(H_R(k))^2;
            den = den + abs(H_L(k))^2;
        end
        
        % full band ILD
        ILD = 10*log10(num/den);
        
        % write data
        if elev < 0
            fprintf(fid,'-%02d %+.3f\n', abs(elev), ILD);
        else
            fprintf(fid,'+%02d %+.3f\n', elev, ILD);
        end
    end
    
    % close file
    fclose(fid);
end
