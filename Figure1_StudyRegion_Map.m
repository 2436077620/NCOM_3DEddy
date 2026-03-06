clc;clear;close all;
% NCOM数据的研究区域：当前是Alaska和Amseas
r = [187-360 308-360 4 65];
r1 = [190-360, 240-360, 36.45, 62.05]; 
r2 = [262-360, 305.1-360, 5, 32.1];
figure;
% 投影类型 'Miller'，'hammer-aitoff','Equidistant Cylindrical'
set(gcf,'outerposition',get(0,'screensize'),'color','w');
set(gca,'fontsize',34,'Fontname','Time new roman','fontweight','bold');hold on;  
m_proj('Miller','longitudes',[r(1) r(2)],'latitudes',[r(3) r(4)]); % 0~360
hold on
% m_map自带水深图
colormap(m_colmap('blues'));  
caxis([-6000 100]);       
m_elev('shadedrelief');
%-----------------------------------海岸线--------------------------------------------  
m_gshhs('ic','patch',[0.3 0.5 0.35],'edgecolor','none');
%-----------------------------------数据边框--------------------------------------------  
m_plot([r1(1) r1(2) r1(2) r1(1) r1(1)], [r1(3) r1(3) r1(4) r1(4) r1(3)],'-','color','r','linewi',4);
m_plot([r2(1) r2(2) r2(2) r2(1) r2(1)], [r2(3) r2(3) r2(4) r2(4) r2(3)],'-','color','r','linewi',4);
m_grid('xtick',7,'ytick',7,'box','on','tickdir','in','LineWidth',1);

% pathname='E:\SWOT_NCOM\';
% print(gcf,'-djpeg','-r450',[pathname,'Figure1.jpg']);