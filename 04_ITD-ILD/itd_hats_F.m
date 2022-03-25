function [ITD] = itd_hats_F(azim,elev)
%-------------------------------------------------------------------------
%   ITD : interaural time difference
%
%   azim : azimuth angle (0 ~ 360 deg) or (-180 ~ +180 deg)
%   elev : elevation angle (-40 ~ +90 deg)
%-------------------------------------------------------------------------
%   Date : July 15, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : retrieval of HATS ITD
%   Synopsis : return HATS ITD
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
root = '../00_Data/ITD_HATS';

% pathname for ITD data
pathname = sprintf('%s/a%03d.txt', root,azim);

% open file
fid = fopen(pathname,'r');
if fid == -1, error('cannot open file : %s',pathname); end

% read data
data = fscanf(fid,'%d %f',[2 inf]);
elevs = data(1,:)';
ITDs = data(2,:)';

% close file
fclose(fid);

% ITD
ITD = ITDs(elevs == elev);
