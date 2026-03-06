clc;clear;close all;
% 绘制delta SSHA 示例图
region   = 'alaska'; % 涡旋所在海域  alaska, amseas
swotpath = 'E:\DATA\SWOTData\'; % SWOT文件存放的路径
% cycle1~7对应num=1；cycle8跨过2023.12.31；cycle9~23对应num=2；cycle24跨过2024.11.19；cycle25~26对应num=3，2025.03.26
cmems_time = '20231231';         % 修改！！！注意！！！（注意下面的n对应的cycle）
slapath  = dir(['E:\DATA\CMEMS_SLAData\','cmems*',cmems_time,'_',region,'.nc']); % 中尺度海面高文件存储路径
outpath  = 'E:\Figure\';
% --------------------------读取CMEMS SLA数据------------------------------
slafile  = [slapath.folder,'\',slapath.name];
lon = double(ncread(slafile,'longitude')); lon(lon<0) = lon(lon<0)+360;
lat = double(ncread(slafile,'latitude'));
[mes.lat,mes.lon] = meshgrid(lat, lon);
mes.time = ncread(slafile,'time')/86400+datenum(1970,1,1);
mes.sla = ncread(slafile,'sla');
clear slafile lon lat
% ---------------------------NCOM模式数据的范围----------------------------
if strcmp(region,'alaska')
    range = [190, 240, 36.45, 62.05]; % 0~360
elseif strcmp(region,'amseas')
    range = [262, 305.1, 5, 32.1];
end
% ----------------------------读取SWOT数据---------------------------------
cycles=dir([swotpath,'cycle*']);
num = 1;
sowtsubfolde = [swotpath,cycles(num).name,'\'];
swotdata = dir([sowtsubfolde,'*.nc']);
i = 1; % nc文件
swotfile=[sowtsubfolde,swotdata(i).name];
swot = ReadSWOT_NCOMrange(swotfile,range); % 读取与模式数据空间范围一致的SWOT，输出结构体变量
SWOTData = struct();
SWOTData.Name = swotdata(i).name;
SWOTData.Time = swot.time;
% 减去中尺度部分（0.125°的网格产品）
[~,mes_day] = min(abs(mes.time - swot.time));
swot.mes_sla = griddata(mes.lon, mes.lat, mes.sla(:,:,mes_day), swot.lon, swot.lat, 'linear');
swot.ssha_sub = swot.ssha - swot.mes_sla;
back = double(isnan(swot.ssha)); back([5,30,35,39,65],:)=1;back(back==0)=nan; 
%% Figure（a）
figure
set(gcf,'outerposition',get(0,'screensize'),'color','w');
set(gca,'fontsize',26,'Fontname','Time new roman','fontweight','bold');hold on; 
m_proj('miller','longitudes',[225.35,227.65],'latitudes',[40.55,42.575]);
hold on
m_pcolor(mes.lon, mes.lat, mes.sla(:,:,mes_day));shading interp;
colormap(m_colmap('diverging'));
c = colorbar;
caxis([0.04 0.16])
c.Label.String = 'SSHA (m)';
m_coast('patch',[0.5 0.5 0.5],'edgecolor',[0.6 0.6 0.6]);
m_grid('xtick',5,'box','fancy','tickdir','out','LineWidth',0.5);
% print(gcf,'-djpeg','-r450',[outpath,'DUACS-SSHA_example1.jpg']);
%% Figure（b）
figure
set(gcf,'outerposition',get(0,'screensize'),'color','w');
set(gca,'fontsize',26,'Fontname','Time new roman','fontweight','bold');hold on; 
m_proj('miller','longitudes',[225.35,227.65],'latitudes',[40.55,42.575]);
hold on
m_pcolor(mes.lon, mes.lat, mes.sla(:,:,mes_day));shading interp;colormap(m_colmap('diverging'));freezeColors;
m_pcolor(swot.lon, swot.lat,back,'FaceAlpha',0.2);colormap(gray);caxis([0 1.5]);freezeColors;
m_pcolor(swot.lon, swot.lat,swot.ssha);shading interp;
colormap(m_colmap('diverging'));
c = colorbar;
caxis([0.04 0.16])
c.Label.String = 'SSHA (m)';
m_coast('patch',[0.5 0.5 0.5],'edgecolor',[0.6 0.6 0.6]);
m_grid('xtick',5,'box','fancy','tickdir','out','LineWidth',0.5);
% print(gcf,'-djpeg','-r450',[outpath,'SWOT-SSHA_example1.jpg']);
%% Figure（c）
load(['E:\SWOT_NCOM\SWOTEddy\',region,'\',swotdata(i).name(1:end-3),'.mat']);
figure;
set(gcf,'outerposition',get(0,'screensize'),'color','w');
set(gca,'fontsize',26,'Fontname','Time new roman','fontweight','bold');hold on; 
m_proj('miller','longitudes',[225.35,227.65],'latitudes',[40.55,42.575]);
hold on
% m_pcolor(mes.lon, mes.lat, mes.sla(:,:,mes_day));shading interp;
m_pcolor(swot.lon, swot.lat,swot.ssha_sub);shading interp;
for i = 1:length(SWOTData.Eddy)
    [lat,lon] = meshgrid(SWOTData.Grid(i).lat,SWOTData.Grid(i).lon);
    if SWOTData.Eddy(i).type == 1
        m_plot(SWOTData.Eddy(i).edge(1,:),SWOTData.Eddy(i).edge(2,:),'-r','linewi',1.5);
    elseif SWOTData.Eddy(i).type == -1 
        m_plot(SWOTData.Eddy(i).edge(1,:),SWOTData.Eddy(i).edge(2,:),'-b','linewi',1.5);
    end
end
colormap(m_colmap('diverging')); 
c = colorbar;
caxis([-0.02 0.03])
c.Label.String = '\Delta SSHA (m)';
m_coast('patch',[0.5 0.5 0.5],'edgecolor',[0.6 0.6 0.6]);
m_grid('xtick',5,'box','fancy','tickdir','out','LineWidth',0.5);
% print(gcf,'-djpeg','-r450',[outpath,'Eddy_example1.jpg']);