function [h2] = ir_window_F(h1)
%-------------------------------------------------------------------------
%   h2 : impulse response (IR) after windowing & zero-padding
%
%   h1 : original IR
%-------------------------------------------------------------------------
%   Date : July 01, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : IR windowing & zero-padding
%   Synopsis : return IR cut out by rectangular window and
%              zero-padded to 512 samples
%	Algorithm : -
%-------------------------------------------------------------------------

% index of start point
i_start = 102;          % 1 ms before the max sample of the ipsilateral IR

% zero-crossing index of IR
[idx] = zci_F(h1);

% index of end point
[~,i_peak] = max(h1);
i_end = i_peak + 120;   % 2.5 ms after the max sample of IR

for n = 1:length(idx)
    if idx(n) >= i_end
        i_end = idx(n);
        break
    end
end

% windowing & zero-padding of IR
h2 = zeros(512,1);
h2(i_start:i_end) = h1(i_start:i_end);
