function [h2] = ir_align_F(h1,idx)
%-------------------------------------------------------------------------
%   h2 : impulse response (IR) after alignment
%
%   h1 : IR before alignment
%   idx : index for alignment
%-------------------------------------------------------------------------
%   Date : June 17, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : alignment of IR at reference point
%   Synopsis : return aligned IR
%	Algorithm : -
%-------------------------------------------------------------------------

% check index
if (idx < 1) || (mod(idx,1) ~= 0)
    error('index must be a natural number');
elseif idx > length(h1)
	error('index must be less than the length of the impulse response');
end

% peak index of h1
[~,i_peak] = max(h1);

if idx == i_peak
    h2 = h1;
elseif idx > i_peak
    h2 = zeros(size(h1));
    di = idx - i_peak;
    h2(di+1:end) = h1(1:end-di);
else
    h2 = zeros(size(h1));
    di = i_peak - idx;
    h2(1:end-di) = h1(di+1:end);
end
