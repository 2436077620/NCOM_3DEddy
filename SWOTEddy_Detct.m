clc;clear;close all;
% 检测SWOT闭合轮廓，用于与NCOM模式数据建立匹配关系
region   = 'amseas'; % 涡旋所在海域  alaska, amseas
swotpath = 'E:\SWOTData\'; % SWOT文件存放的路径
% cycle1~7对应num=1；cycle8跨过2023.12.31；cycle9~23对应num=2；cycle24跨过2024.11.19；cycle25~26对应num=3，2025.03.26
cmems_time = '20250326';         % 修改！！！注意！！！（注意下面的n对应的cycle）
slapath  = dir(['E:\CMEMS_SLAData\','cmems*',cmems_time,'_',region,'.nc']); % 中尺度海面高文件存储路径
outpath  = ['E:\SWOT_NCOM\SWOTEddy\',region,'\'];
% ---------------------------NCOM模式数据的范围----------------------------
if strcmp(region,'alaska')
    range = [190, 240, 36.45, 62.05]; % 0~360
elseif strcmp(region,'amseas')
    range = [262, 305.1, 5, 32.1];
end
% --------------------------读取CMEMS SLA数据------------------------------
slafile  = [slapath.folder,'\',slapath.name];
lon = double(ncread(slafile,'longitude')); lon(lon<0) = lon(lon<0)+360;
lat = double(ncread(slafile,'latitude'));
[mes.lat,mes.lon] = meshgrid(lat, lon);
mes.time = ncread(slafile,'time')/86400+datenum(1970,1,1);
mes.sla = ncread(slafile,'sla');
clear slafile lon lat
%% ----------------------------读取SWOT数据---------------------------------
cycles=dir([swotpath,'cycle*']);
for n = 30:30 % cycle：1~7半年的数据
% n = 1; 
    sowtsubfolde = [swotpath,cycles(n).name,'\'];
    swotdata = dir([sowtsubfolde,'*.nc']);
    for i =1:257%length(swotdata) % cycle8分为1:313 & 314:length(swotdata); cycle24分为1:203 & 204:length(swotdata); cycle30分为1:257 & 258:length(swotdata)
%     i = 3; % nc文件
        swotfile=[sowtsubfolde,swotdata(i).name];
        swot = ReadSWOT_NCOMrange(swotfile,range); % 读取与模式数据空间范围一致的SWOT，输出结构体变量
        SWOTData = struct();
        if isfield(swot,'ssha') % 裁剪后存在ssha变量，即经过指定范围，进行下一步
            if mean(isnan(swot.ssha),'all')==1 % ssha全是nan，则跳过
                continue;
            end
            SWOTData.Name = swotdata(i).name;
            SWOTData.Time = swot.time;
% ----------------------------检测SWOT涡旋---------------------------------
            % 减去中尺度部分（0.125°的网格产品）
            [~,mes_day] = min(abs(mes.time - swot.time));
            swot.mes_sla = griddata(mes.lon, mes.lat, mes.sla(:,:,mes_day), swot.lon, swot.lat, 'linear');
            swot.ssha_sub = swot.ssha - swot.mes_sla;
            % 1 选择强振幅；2 直径能够被NCOM数据包括的涡旋
            SWOTData.Eddy = SWOT_Eddydec(swot.lon,swot.lat,swot.ssha_sub,0.4,8); %振幅>0.4cm，直径>8km (耗时长)
            if isfield(SWOTData.Eddy,'center') % 如果检测出有涡旋，进行下一步
% ------------------裁剪检测出的SWOT Eddy的海面高网格----------------------
                for k =1:length(SWOTData.Eddy)
                    SWOTData.Grid(k) = cutswotEddy(swot,SWOTData.Eddy(k)); % 两倍空间范围的网格，裁剪并插值到规则网格
                end
                save([outpath,SWOTData.Name(1:end-3),'.mat'],'SWOTData');
            end
        end
        disp(['cycle',num2str(n),': ',num2str(i),'/',num2str(length(swotdata))])
    end
end
%% 
%{
figure
m_proj('miller','longitudes',[188,242],'latitudes',[33,65]);
hold on
% m_pcolor(mes.lon, mes.lat, mes.sla(:,:,mes_day));shading interp;
m_pcolor(swot.lon, swot.lat,swot.ssha);shading interp;
colormap(m_colmap('diverging')); 
c.Label.String = 'SSHA (m)';
m_coast('patch',[0.5 0.5 0.5],'edgecolor',[0.6 0.6 0.6]);
m_grid('xtick',5,'box','fancy','tickdir','out','fonts',16,'LineWidth',0.5);
% 绘图提取之后的涡旋
clc;clear;
load('E:\SWOT_NCOM\SWOTEddy\alaska\SWOT_L3_LR_SSH_Expert_001_166_20230727T030232_20230727T035359_v1.0.mat');
figure;
m_proj('miller','longitudes',[188,200],'latitudes',[35,58]);
hold on
% m_pcolor(mes.lon, mes.lat, mes.sla(:,:,mes_day));shading interp;
for i = 1:length(SWOTData.Eddy)
    [lat,lon] = meshgrid(SWOTData.Grid(i).lat,SWOTData.Grid(i).lon);
    m_pcolor(lon,lat,SWOTData.Grid(i).ssha_sub);shading interp;
    if SWOTData.Eddy(i).type == 1
        m_plot(SWOTData.Eddy(i).edge(1,:),SWOTData.Eddy(i).edge(2,:),'-r');
    elseif SWOTData.Eddy(i).type == -1 
        m_plot(SWOTData.Eddy(i).edge(1,:),SWOTData.Eddy(i).edge(2,:),'-b');
    end
end
colormap(m_colmap('diverging')); 
c.Label.String = 'SSHA (m)';
m_coast('patch',[0.5 0.5 0.5],'edgecolor',[0.6 0.6 0.6]);
m_grid('xtick',5,'box','fancy','tickdir','out','fonts',16,'LineWidth',0.5);
%}