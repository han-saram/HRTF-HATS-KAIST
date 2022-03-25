function [h_L,h_R] = hrir_hats_F(azim,elev,mode)
%-------------------------------------------------------------------------
%   h_L : Left HRIR of HATS
%   h_R : Right HRIR of HATS
%
%   azim : azimuth angle (0 ~ 360 deg) or (-180 ~ +180 deg)
%   elev : elevation angle (-40 ~ +90 deg)
%   mode : 1 (1 ms delayed HRIR) or any-number (raw HRIR)
%-------------------------------------------------------------------------
%   Date : July 05, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : retrieval of HATS HRIR
%   Synopsis : return stereo HRIR of HATS
%	Algorithm : -
%-------------------------------------------------------------------------

if nargin == 2, mode = 1; end

% check azimuth
if mod(azim,5) ~= 0
    error('azimuth must be a multiple of 5');
elseif (azim < -180) || (azim > 360)
    error('azimuth must be 0 ~ 360 deg or -180 ~ +180 deg');
elseif (azim >= -180) && (azim < 0)
    azim = 360 + azim;
elseif azim == 360
    azim = 0;
end

% check elevation
if mod(elev,5) ~= 0
    error('elevation must be a multiple of 5');
elseif (elev < -40) || (elev > 90)
	error('elevation must be -40 ~ +90 deg');
end

% root of directory
if mode == 1
    root = '../00_Data/HRTF_HATS';
else
    root = '../00_Data/HRTF_HATS_Raw';
end

% pathname for IR data
if elev >= 0
    pathname = sprintf('%s/a%03d/a%03de+%02d.txt', root,azim,azim,elev);
else
    pathname = sprintf('%s/a%03d/a%03de-%02d.txt', root,azim,azim,abs(elev));
end

% open file
fid = fopen(pathname,'r');
if fid == -1, error('cannot open file : %s',pathname); end

% read data
h = fscanf(fid,'%f %f',[2 inf]);
h_L = h(1,:)';
h_R = h(2,:)';

% close file
fclose(fid);
