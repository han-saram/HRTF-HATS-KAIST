function [h,Fs] = ir_hats_F(azim,elev,root)
%-------------------------------------------------------------------------
%   h : stereo impulse responses (IR) of HATS
%   Fs : sampling frequency
%
%   azim : azimuth angle (0 ~ 360 deg) or (-180 ~ +180 deg)
%   elev : elevation angle (-40 ~ +90 deg)
%   root : root of directory
%-------------------------------------------------------------------------
%   Date : June 03, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : pair formation of HATS IR
%   Synopsis : return stereo IR of HATS in first two columns of h
%              such that left is first column, right is second column
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

% pathname for IR data
if elev >= 0
    pathname_L = sprintf('%s/a%03d/La%03de+%02d.txt', root,azim,azim,elev);
    pathname_R = sprintf('%s/a%03d/Ra%03de+%02d.txt', root,azim,azim,elev);
else
    pathname_L = sprintf('%s/a%03d/La%03de-%02d.txt', root,azim,azim,abs(elev));
    pathname_R = sprintf('%s/a%03d/Ra%03de-%02d.txt', root,azim,azim,abs(elev));
end

% open file
fid_L = fopen(pathname_L,'r');
fid_R = fopen(pathname_R,'r');
if fid_L == -1
	error('cannot open file : %s',pathname_L);
elseif fid_R == -1
    error('cannot open file : %s',pathname_R);
end
fgetl(fid_L);
fgetl(fid_R);

% read data
h_L = fscanf(fid_L,'%f %f %f',[3 inf]);
h_R = fscanf(fid_R,'%f %f %f',[3 inf]);
h(:,1) = h_L(2,:)';
h(:,2) = h_R(2,:)';
Fs = 1/h_L(1,2);

% close file
fclose(fid_L);
fclose(fid_R);
