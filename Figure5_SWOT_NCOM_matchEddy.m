clc;clear;close all
% SWOT Eddy和NCOM Eddy时空匹配的图像 (数据：Eddy_OW_Dataset)
% 用于Model validation
region    = 'amseas'; % 涡旋所在海域  alaska, amseas
eddy_path = 'E:\DATA\Eddy_OW_Dataset\';
eddy_file = dir([eddy_path,'*.mat']);
n=1;
load([eddy_path,eddy_file(n).name]);
[Lat1, Lon1] = meshgrid(Input_grid.lat,Input_grid.lon);
[Lat2, Lon2] = meshgrid(Output_grid.lat,Output_grid.lon);
disp(eddy_file(n).name);
center= cell2mat({Output_eddy.center}');% NCOM涡心
dis = distance(center(:,2),center(:,1),Input_eddy.center(2),Input_eddy.center(1))./180.*pi.*6371;% 与SWOT涡心距离
%% SWOT涡旋对应的完整数据
range = [262, 305.1, 5, 32.1];
% ------读取CMEMS SLA数据-------
cmems_time = '20250326';         % 修改！！！注意！！！（注意下面的n对应的cycle）
slapath  = dir(['E:\DATA\CMEMS_SLAData\','cmems*',cmems_time,'_',region,'.nc']); % 中尺度海面高文件存储路径
slafile  = [slapath.folder,'\',slapath.name];
lon = double(ncread(slafile,'longitude')); lon(lon<0) = lon(lon<0)+360;
lat = double(ncread(slafile,'latitude'));
[mes.lat,mes.lon] = meshgrid(lat, lon);
mes.time = ncread(slafile,'time')/86400+datenum(1970,1,1);
mes.sla = ncread(slafile,'sla');
clear slafile lon lat
% ------------SWOT-------------
nm = eddy_file(n).name;
swotfile = dir(['E:\DATA\SWOTData\cycle_',nm(6:8),'\*',nm(5:13),'*.nc']);
swot = ReadSWOT_NCOMrange([swotfile.folder,'\',swotfile.name],range); % 读取与模式数据空间范围一致的SWOT，输出结构体变量
load(['E:\DATA\SWOTEddy\',swotfile.name(1:end-3),'.mat']);
[~,mes_day] = min(abs(mes.time - swot.time));% 减去中尺度部分（0.125°的网格产品）
swot.mes_sla = griddata(mes.lon, mes.lat, mes.sla(:,:,mes_day), swot.lon, swot.lat, 'linear');
swot.ssha_sub = swot.ssha - swot.mes_sla;
back = ones(size(swot.ssha));
%% -----------画图-------------
figure;
set(gcf,'position',[200,50,1100,700],'color','w');
set(gca,'fontsize',26,'Fontname','Time new roman','fontweight','bold');hold on; 
m_proj('miller','longitudes',[284.8,287],'latitudes',[20.81,22.45]);
hold on
m_pcolor(mes.lon, mes.lat, inpaint_nans(mes.sla(:,:,mes_day)).*100);shading interp;colormap(m_colmap('diverging'));caxis([-3 13]);freezeColors;
m_pcolor(swot.lon, swot.lat,back,'FaceAlpha',0.2);shading interp;colormap(gray);caxis([0 1.5]);freezeColors;
m_pcolor(swot.lon, swot.lat,swot.ssha.*100);shading interp;
colormap(m_colmap('diverging')); 
c = colorbar;
caxis([-3 8])
c.Label.String = 'SSHA (cm)';
m_plot([Input_grid.lon(1) Input_grid.lon(1) Input_grid.lon(end) Input_grid.lon(end) Input_grid.lon(1)],...
    [Input_grid.lat(2) Input_grid.lat(end-1) Input_grid.lat(end-1) Input_grid.lat(2) Input_grid.lat(2)],'-k','linewi',1.5)
m_gshhs_i('patch',[0.5 0.5 0.5],'edgecolor',[0.6 0.6 0.6]);
m_grid('xtick',5,'box','fancy','tickdir','out','LineWidth',0.5);
% print(gcf,'-djpeg','-r450',['E:\Figure\All_',eddy_file(n).name(1:end-4),'.jpg']);
%% 输入网格的图像
[X,Y] = LL3XY(Lon1, Lat1);
[X0,Y0] = LL3XY(Input_eddy.center(1), Input_eddy.center(2));
[Edge(1,:),Edge(2,:)] = LL3XY(Input_eddy.edge(1,:),Input_eddy.edge(2,:));
% 涡旋轮廓
edgelen = length(Edge(1,:));
edgeLon = interp1(1:edgelen,Edge(1,:),1:0.05:edgelen,'pichp');
edgeLat = interp1(1:edgelen,Edge(2,:),1:0.05:edgelen,'pichp');
figure
set(gcf,'position',[200,50,620,500],'color','w');
set(gca,'position',[0.15,0.17,0.63,0.75],'FontName','Times New Roman','FontSize',22,'FontWeight','bold');hold on;
pcolor(X-X0, Y-Y0, inpaint_nans(Input_grid.ssha_sub,1).*100);shading interp;
colormap(flip(othercolor('RdBu11')));
c = colorbar('location','east','position',[0.8 0.17 0.03 0.75]);
c.Label.String = '\Delta SSHA (cm)';
caxis([2.1 3.7]);
plot(edgeLon-X0,edgeLat-Y0,'-','color',[.4 .4 .4],'linewi',1.5);
xlim([-12 12]); ylim([-12 12]); xlabel('dX (km)'); ylabel('dY (km)');
plot([-12 -12 12 12 -12],[-12 12 12 -12 -12],'-k','linewi',2)
grid on
% print(gcf,'-djpeg','-r450',['E:\Figure\',eddy_file(n).name(1:end-4),'.jpg']);
%% 输出网格的水下三维
posit = [700, 0.63, 0.82]; % gcf宽度，gca右边界，colorbar左边界
var_name = 'vorticity';
var_str = '\zeta/f';
var_range= [-0.4,0.4];
figure;%-------------------------------------------------
set(gcf,'position',[200,50,posit(1),950],'color','w');
set(gca,'position',[0.17,0.08,posit(2),0.9],'FontName','Times New Roman','FontSize',22,'FontWeight','bold');hold on;
lonlim = []; latlim = []; eddycen = [];
%--------------绘制参数-----------------
for dd=1:length(Output_eddy)-2
    eddy_para = Output_eddy(dd); % 被追踪的涡旋
    dep = -eddy_para.depth; %深度
    [X,Y] = LL3XY(eddy_para.center(1),eddy_para.center(2));
    eddycen(dd,:) = [X, Y, dep];
    if ismember(dd,2) % 仅绘制特定层的图像
        continue;
    end
    d0 = find(Output_grid.depth==-dep);
    lonlim = [min([eddy_para.edgelon, lonlim]),max([eddy_para.edgelon, lonlim])]; % 范围
    latlim = [min([eddy_para.edgelat, latlim]),max([eddy_para.edgelat, latlim])]; % 范围
    [lat_f,lon_f] = meshgrid(Output_grid.lat,Output_grid.lon);
    var0 = Output_grid.vorticity(:,:,d0)./gsw_f(lat_f); % 绘制的参数（相对涡度需要基于f归一化）
    [xlen,ylen] = size(var0);
    % ----插值----
    % 层级参数
    res =0.02; % 提升分辨率
    lon = interp1(1:xlen,Output_grid.lon,1:res:xlen);
    lat = interp1(1:ylen,Output_grid.lat,1:res:ylen);
    [lat_re,lon_re] = meshgrid(lat,lon);
    var = interp2(lat_f,lon_f,var0,lat_re,lon_re,'linear');
    clear xlen ylen lon lat var0
    % 涡旋轮廓
    res2 = 0.1;
    edgelen = length(eddy_para.edgelon);
    edgelon = interp1(1:edgelen,eddy_para.edgelon,1:res2:edgelen,'pichp');
    edgelat = interp1(1:edgelen,eddy_para.edgelat,1:res2:edgelen,'pichp');
    % ---涡旋内---
    in=inpolygon(lon_re, lat_re, edgelon, edgelat); %内部点
    % ----画图----
    [x1,y1] = LL3XY(lon_re,lat_re);
    [x2,y2] = LL3XY(lon_f, lat_f);
    surf(x1-eddycen(1,1),y1-eddycen(1,2),dep.*ones(size(lat_re)),var); shading interp% 在z=h的位置绘制图层
    [e1,e2] = LL3XY(edgelon,edgelat);
    plot3(e1-eddycen(1,1), e2-eddycen(1,2), repmat(dep,1,(edgelen-1)/res2+1),'-k','linewi',0.6,'color',[.4 .4 .4]);
    clear dep lat_re lon_re var edgelen edgelon edgelat in lat_f lon_f in0
end
colormap(m_colmap('diverging'));
h1=colorbar('location','east','position',[posit(3) 0.19 0.025 0.72]);
h1.Label.String = var_str;
caxis(var_range);
% 涡旋中心的垂向连线
eddycen0 = [eddycen(:,1)-eddycen(1,1), eddycen(:,2)-eddycen(1,2)];
plot3(eddycen0(:,1), eddycen0(:,2), eddycen(:,3), '--' ,'linewi',1.5, 'color', [102 109 101]./255);
plot3(eddycen0(:,1), eddycen0(:,2), eddycen(:,3)+0.2, '.','MarkerSize',20,'color', [102 109 101]./255)
view(-60,10);
set(gca,'xlim',[-10,14],'ylim',[-10,11],'zlim',[eddycen(end,3)-0.1,0.2],'color','w')
pos = axis;
xlabel('dX (km)','rotation',40, 'position',[pos(1)+0.08*(pos(2)-pos(1)), pos(3)-0.3*(pos(4)-pos(3)), pos(5)]);
ylabel('dY (km)', 'rotation',-14, 'position',[pos(1)-0.32*(pos(2)-pos(1)), pos(3)+0.22*(pos(4)-pos(3)), pos(5)]);
zlabel('Depth (m)');
box on
grid on
% print(gcf,'-djpeg','-r450',['E:\Figure\NCOM_',var_name,'_swot_025_007.jpg']);