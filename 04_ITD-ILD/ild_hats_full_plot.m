%-------------------------------------------------------------------------
%   Date : July 16, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : HATS full band ILD profile & map
%   Synopsis : plot HATS full band ILD profile & map
%	Algorithm : -
%-------------------------------------------------------------------------

clc
clear
close all

% azimuth angles (-180 ~ +180 deg)
azims = (-180:5:180)';
N_azim = length(azims);

% elevation angles (-40 ~ +90 deg)
elevs = (-40:5:90)';
N_elev = length(elevs);

% initialization
ILD_map = zeros(N_elev,N_azim);

% ILD Map
for j = 1:N_azim
    azim = azims(j);
    
    for i = 1:N_elev
        elev = elevs(i);
        
        % ILD retrieval
        [ILD] = ild_hats_full_F(azim,elev);
        ILD_map(i,j) = ILD;
    end
end

% plot
figure
plot(azims,ILD_map(elevs == 0,:),'LineWidth',1.2);
hold on
plot(azims,ILD_map(elevs == 15,:),'LineWidth',1.2);
hold on
plot(azims,ILD_map(elevs == 30,:),'LineWidth',1.2);
hold on
plot(azims,ILD_map(elevs == 45,:),'LineWidth',1.2);
hold on
plot(azims,ILD_map(elevs == 60,:),'LineWidth',1.2);
hold on
plot(azims,ILD_map(elevs == 75,:),'LineWidth',1.2);
hold on
plot(azims,ILD_map(elevs == 90,:),'LineWidth',1.2);
legend('Elev. 0\circ','Elev. 15\circ','Elev. 30\circ','Elev. 45\circ','Elev. 60\circ','Elev. 75\circ','Elev. 90\circ','Location','northwest');
xlabel('Azimuth [deg.]','fontsize',12); ylabel('Level [dB]','fontsize',12);
title('ILD profile of HATS - full band','fontsize',14);
axis([azims(1) azims(end) -20 20]);
grid on
set(gcf,'position',[100 500 700 500]);

figure
imagesc(azims,elevs,ILD_map); axis xy
colormap jet; colorbar
xlabel('Azimuth [deg.]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
title('ILD map of HATS - full band','fontsize',14);
title(colorbar,'Level [dB]')
axis([azims(1) azims(end) elevs(1) elevs(end)]);
grid on
set(gcf,'position',[900 500 700 500]);
