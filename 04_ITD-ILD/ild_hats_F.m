function [ILD,f] = ild_hats_F(azim,elev)
%-------------------------------------------------------------------------
%   ILD : interaural level difference
%   f : freq. axis
%
%   azim : azimuth angle (0 ~ 360 deg) or (-180 ~ +180 deg)
%   elev : elevation angle (-40 ~ +90 deg)
%-------------------------------------------------------------------------
%   Date : July 16, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : retrieval of HATS ILD
%   Synopsis : return HATS ILD
%	Algorithm : -
%-------------------------------------------------------------------------

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
root = '../00_Data/ILD_HATS';

% pathname for ILD data
if elev < 0
    pathname = sprintf('%s/a%03d/a%03de-%02d.txt', root,azim,azim,abs(elev));
else
    pathname = sprintf('%s/a%03d/a%03de+%02d.txt', root,azim,azim,elev);
end

% open file
fid = fopen(pathname,'r');
if fid == -1, error('cannot open file : %s',pathname); end

% read data
data = fscanf(fid,'%d %f',[2 inf]);
f = data(1,:)';
ILD = data(2,:)';

% close file
fclose(fid);
