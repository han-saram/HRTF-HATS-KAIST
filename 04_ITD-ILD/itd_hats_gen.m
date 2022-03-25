%-------------------------------------------------------------------------
%   Date : July 15, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : HATS HRTF ITD calculation
%   Synopsis : calculate ITDs from HATS HRTF and store them in text files
%	Algorithm : -
%-------------------------------------------------------------------------

clc
clear
close all

% make path
path = '../00_Data/ITD_HATS';
mkdir(path);

% minimum phase lowpass filter design (fc = 1.5 kHz)
Fs = 48e3;
Fpass = 1.0e3;
Fstop = 2.5e3;
Apass = 0.01;
Astop = 80;
LPF_spec = fdesign.lowpass('Fp,Fst,Ap,Ast',Fpass,Fstop,Apass,Astop,Fs);
LPF = design(LPF_spec,'equiripple','minphase',true,'SystemObject',true);
%fvt = fvtool(LPF,'Fs',Fs,'Color','White');

% azimuth angles (0 ~ 355 deg)
for azim = 0:5:355
    fprintf('calculating ITDs at azimuth %d deg\n',azim);
    
    % pathname
    name = sprintf('a%03d.txt', azim);
    pathname = [path '/' name];
    
    % open file
    fid = fopen(pathname,'w');
    
    % elevation angles (-40 ~ +90 deg)
    for elev = -40:5:90
        % HRIR retrieval
        [h_L0,h_R0] = hrir_hats_F(azim,elev);
        
        % lowpass filtering
        h_L1 = LPF(h_L0);
        h_R1 = LPF(h_R0);
        
        % upsampling: 48 kHz x 4 = 192 kHz
        p = 4; q = 1;
        Fs2 = Fs*p/q;
        h_L = resample(h_L1,p,q);
        h_R = resample(h_R1,p,q);
        
        % interaural cross-correlation
        tau_max = 1000e-6;
        lag_max = tau_max*Fs2;
        [R_lr,lags] = xcorr(h_L,h_R,lag_max,'normalized');
        tau = lags'/Fs2;
        
        % ITD
        [~,idx] = max(R_lr);
        ITD = tau(idx);
        
        % write data
        if elev < 0
            fprintf(fid,'-%02d %+.12f\n', abs(elev), ITD);
        else
            fprintf(fid,'+%02d %+.12f\n', elev, ITD);
        end
    end
    
    % close file
    fclose(fid);
end
