%-------------------------------------------------------------------------
%   Date : June 30, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------

clc
clear
close all

% root of directory
root = '../00_Data/TF_HATS';

% load IR
[L_h,Fs] = ir_hats_F(-90,0,root);
[R_h,~] = ir_hats_F(90,0,root);
L_ips = L_h(:,1);
R_ips = R_h(:,2);

% sample numbers of max amplitude
[m_L,i_L] = max(L_ips);
[m_R,i_R] = max(R_ips);

% check the sample numbers
if i_L == i_R
    i_M = i_L;
    t_M = 1000*(i_M-1)/Fs;      % ms
else
	error('sample numbers of max amplitude are different');
end

% sample number of start point
i_sp = i_M - 48;                % 1 ms before the max amplitude
t_sp = 1000*(i_sp-1)/Fs;        % ms

% time axis
N = length(R_ips);
time_axis = 1000*(0:N-1)/Fs;    % ms

% plot
figure
plot(time_axis,L_ips);
hold on
plot(time_axis,R_ips);
hold on
xline(t_M,'-',{[num2str(i_M),'th sample'],[num2str(t_M),' ms']});
hold on
xline(t_sp,'-',{[num2str(i_sp),'nd sample'],[num2str(t_sp),' ms']});
legend('Left ipsilateral','Right ipsilateral');
xlabel('Time [ms]','fontsize',12); ylabel('Amplitude','fontsize',12);
title('Impulse Responses of Ipsilateral Directions','fontsize',14);
axis([0 10 -0.4 0.4]);
grid on
set(gcf,'position',[300 500 900 350]);
