function [h,Fs] = ir_ref1_F(elev,mic_dir,root)
%-------------------------------------------------------------------------
%   h : impulse response (IR) at reference point
%   Fs : sampling frequency
%
%   elev : elevation angle (-40 ~ +90 deg)
%   mic_dir : mic direction [-40, -30, 0, +30, +60, +90 deg]
%   root : root of directory
%-------------------------------------------------------------------------
%   Date : June 17, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : formation of IR at reference point
%   Synopsis : return IR at ref. point measured by Mic-R in 1 direction
%	Algorithm : -
%-------------------------------------------------------------------------

% check elevation
if mod(elev,5) ~= 0
    error('elevation must be a multiple of 5');
elseif (elev < -40) || (elev > 90)
	error('elevation must be -40 ~ +90 deg');
end

% check mic direction
if (mic_dir~=-40)&&(mic_dir~=-30)&&(mic_dir~=0)&&(mic_dir~=30)&&(mic_dir~=60)&&(mic_dir~=90)
    error('mic direction must be one of -40, -30, 0, +30, +60, +90 deg');
end

% pathname for IR data
if mic_dir < 0
    if elev < 0
        pathname = sprintf('%s/Mic-R-%02d/Ref_e-%02d [Mic-R-%02d SPK-%02d].txt', root,abs(mic_dir),abs(elev),abs(mic_dir),abs(elev));
    else
        pathname = sprintf('%s/Mic-R-%02d/Ref_e+%02d [Mic-R-%02d SPK+%02d].txt', root,abs(mic_dir),elev,abs(mic_dir),elev);
    end
else
    if elev < 0
        pathname = sprintf('%s/Mic-R+%02d/Ref_e-%02d [Mic-R+%02d SPK-%02d].txt', root,mic_dir,abs(elev),mic_dir,abs(elev));
    else
        pathname = sprintf('%s/Mic-R+%02d/Ref_e+%02d [Mic-R+%02d SPK+%02d].txt', root,mic_dir,elev,mic_dir,elev);
    end
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
