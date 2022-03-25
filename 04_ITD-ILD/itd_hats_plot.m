%-------------------------------------------------------------------------
%   Date : July 15, 2021
%   Place : Human Lab, KAIST
%   Coder : Lee, Gyeong Tae
%-------------------------------------------------------------------------
%	Title : HATS ITD profile & map
%   Synopsis : plot HATS ITD profile & map
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
ITD_map = zeros(N_elev,N_azim);

% ITD Map
for j = 1:N_azim
    azim = azims(j);
    
    for i = 1:N_elev
        elev = elevs(i);
        
        % ITD retrieval
        [ITD] = itd_hats_F(azim,elev);
        ITD_map(i,j) = ITD;
    end
end
ITD_map = ITD_map*1e6;      % micro-seconds

% plot
figure
plot(azims,ITD_map(elevs == 0,:),'LineWidth',1.2);
hold on
plot(azims,ITD_map(elevs == 15,:),'LineWidth',1.2);
hold on
plot(azims,ITD_map(elevs == 30,:),'LineWidth',1.2);
hold on
plot(azims,ITD_map(elevs == 45,:),'LineWidth',1.2);
hold on
plot(azims,ITD_map(elevs == 60,:),'LineWidth',1.2);
hold on
plot(azims,ITD_map(elevs == 75,:),'LineWidth',1.2);
hold on
plot(azims,ITD_map(elevs == 90,:),'LineWidth',1.2);
legend('Elev. 0\circ','Elev. 15\circ','Elev. 30\circ','Elev. 45\circ','Elev. 60\circ','Elev. 75\circ','Elev. 90\circ','Location','northwest');
xlabel('Azimuth [deg.]','fontsize',12); ylabel('Time [\mus]','fontsize',12);
title('ITD profile of HATS','fontsize',14);
axis([azims(1) azims(end) -800 800]);
grid on
set(gcf,'position',[100 500 700 500]);

figure
imagesc(azims,elevs,ITD_map); axis xy
colormap jet;
xlabel('Azimuth [deg.]','fontsize',12); ylabel('Elevation [deg.]','fontsize',12);
title('ITD map of HATS','fontsize',14);
title(colorbar,'Time [\mus]')
axis([azims(1) azims(end) elevs(1) elevs(end)]);
grid on
set(gcf,'position',[900 500 700 500]);
