function [h,Fs] = ir_ref5_F(elev,root)
%-------------------------------------------------------------------------
%   h : impulse response (IR) at reference point
%   Fs : sampling frequency
%
%   elev : elevation angle (-40 ~ +90 deg)
%   root : root of directory
%-------------------------------------------------------------------------
%   Date : June 17, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : formation of IR at reference point
%   Synopsis : return IR at ref. point measured by Mic-R in 5 directions
%	Algorithm : -
%-------------------------------------------------------------------------

% check elevation
if mod(elev,5) ~= 0
    error('elevation must be a multiple of 5');
elseif (elev < -40) || (elev > 90)
	error('elevation must be -40 ~ +90 deg');
end

% pathname for IR data
if (elev >= -40) && (elev <= -15)
    pathname = sprintf('%s/Ref_e-%02d [Mic-R-30 SPK-%02d].txt', root,abs(elev),abs(elev));
elseif (elev >= -10) && (elev <= -5)
    pathname = sprintf('%s/Ref_e-%02d [Mic-R+00 SPK-%02d].txt', root,abs(elev),abs(elev));
elseif (elev >= 0) && (elev <= 15)
    pathname = sprintf('%s/Ref_e+%02d [Mic-R+00 SPK+%02d].txt', root,elev,elev);
elseif (elev >= 20) && (elev <= 45)
    pathname = sprintf('%s/Ref_e+%02d [Mic-R+30 SPK+%02d].txt', root,elev,elev);
elseif (elev >= 50) && (elev <= 75)
    pathname = sprintf('%s/Ref_e+%02d [Mic-R+60 SPK+%02d].txt', root,elev,elev);
else
    pathname = sprintf('%s/Ref_e+%02d [Mic-R+90 SPK+%02d].txt', root,elev,elev);
end

% open file
fid = fopen(pathname,'r');
if fid == -1
	error('cannot open file : %s',pathname);
end
fgetl(fid);

% read data
h_raw = fscanf(fid,'%f %f %f',[3 inf]);
h = h_raw(2,:)';
Fs = 1/h_raw(1,2);

% close file
fclose(fid);
