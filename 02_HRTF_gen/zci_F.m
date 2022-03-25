function [idx] = zci_F(x,mode)
%-------------------------------------------------------------------------
%   idx : index for which x crosses zero
%
%   x : time signal
%   mode : 1(upward) or -1(downward) or any-number(all zero-crossing)
%-------------------------------------------------------------------------
%   Date : July 01, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : zero-crossing index
%   Synopsis : return indices of zero-crossing points
%	Algorithm : -
%-------------------------------------------------------------------------

if nargin == 1, mode = 0; end

x = x(:);
N = length(x);
x1 = x(1:N-1);
x2 = x(2:N);

switch mode
    case 1
        idx = find(x1<=0 & x2>0);
    case -1
        idx = find(x1>=0 & x2<0);
    otherwise
        idx = [find(x1<=0 & x2>0); find(x1>=0 & x2<0)];
        idx = sort(idx);
end

if isempty(idx), error('there is no zero-crossing point'); end
