clc;clear;close all;
% 选择样本，读取并绘制从NCOM数据中识别涡旋的过程（以海表面2D场为例）
% 方法：OW参数法 1）涡度计算 2）OW参数计算 3）轮廓提取 4）筛选涡旋
% 包括amseas全局海域、大样本海域、（中尺度涡及其周围）小样本海域
model_path= 'E:\DATA\NCOM_model\';     % 储存模式数据的路径
eddypath = 'E:\DATA\NCOMGrid\';
res = 3.3;
hpkm = 50;  % 对模式数据高通滤波的尺寸
eddy_file = '_001_173_';
% ---------------------------NCOM模式数据的范围----------------------------
range = [262, 305.1, 5, 32.1];
figrange = [range(1:2);range(3:4)];
eddydata = dir([eddypath,'*',eddy_file,'*.mat']);
load([eddypath,eddydata.name]);
Model = ReadNCOMmodel([model_path, ModelData.NCOM_Name],res,hpkm);% 读model文件
%% NCOM Vorticity & Okubo-Weiss 参数计算
[lat_m, lon_m] = meshgrid(Model.lat, Model.lon);
surf_el_hp = filt2(Model.surf_el,res,hpkm,'hp');
    layer = 1;
    u_hp = filt2(Model.water_u(:,:,layer),res,hpkm,'hp');
    v_hp = filt2(Model.water_v(:,:,layer),res,hpkm,'hp');
    % --------------------相对涡度----------------------------
    [dudx,dudy] = gradient(permute(u_hp,[2,1,3]),(1/30)/180*pi*6371*10^3); % u
    [dvdx,dvdy] = gradient(permute(v_hp,[2,1,3]),(1/30)/180*pi*6371*10^3); % v
    vorticity = permute(dvdx-dudy,[2,1,3]); % 逆时针为+，顺时针为-
    zeta = vorticity./gsw_f(lat_m);% 利用f归一化   
    % ----------------------应变率----------------------------
    Sn = permute(dudx-dvdy,[2,1,3]); %应变的法向分量
    Ss = permute(dvdx+dudy,[2,1,3]); %应变的剪切分量
    W = Sn.^2 + Ss.^2 - vorticity.^2;% 计算OW参数
    W_norm = zscore_normalize(W); %(基于全场的标准差)归一化涡度的范围选取！！！
    clear dudx dudy dvdx dvdy
    W_noiseless = filt2(W_norm, 1/30, 1/12, 'lp');
    OWValue = -0.2;
    ShapeError = 0.55;
    %% ----------------OW参数---------------------------
    figure%('visible','off')
    set(gcf,'outerposition',get(0,'screensize'),'color','w');
    set(gca,'fontsize',26,'Fontname','Time new roman','fontweight','bold');hold on; 
    m_proj('miller','longitudes',[360-72 360-66],'latitudes',[12,18]);
    hold on
    m_pcolor(lon_m, lat_m, W_norm);shading interp;
    colormap(flip(othercolor('RdBu11')));
    c = colorbar('fontsize',26);
    c.Label.String = 'Okubo–Weiss parameter';
    caxis([-0.601 0.6])
    m_coast('patch',[0.5 0.5 0.5],'edgecolor',[0.2 0.2 0.2]);
    m_grid('xtick',5,'box','fancy','tickdir','out','LineWidth',0.5);
%     print(gcf,'-djpeg','-r450',['E:\Figure\','OW_surf_raw_all.jpg']);
    %% 平滑滤波和提取轮廓（阈值0W<-0.2）
    figure('visible','off');
    [c,~]=contour(lon_m, lat_m, W_noiseless, [OWValue 9999]);
    s = getcontourlines(c); % getcontourlines 解析contour得到的"闭合"轮廓
    ss=0;
    k = struct();
    warning off
    for i = 1:length(s) 
    % 保留内部数值<-0.2的轮廓线(Williams et al. 2011)
        in = inpolygon(lon_m, lat_m, s(i).x, s(i).y); %内部点
        if mean(W_noiseless(in),'omitnan') <= OWValue  %内部值<-0.2
            ss=ss+1;
            k(ss).edgelon = s(i).x;
            k(ss).edgelat = s(i).y;
        end
    end
    figure%('visible','off')
    set(gcf,'outerposition',get(0,'screensize'),'color','w');
    set(gca,'fontsize',26,'Fontname','Time new roman','fontweight','bold');hold on; 
    m_proj('miller','longitudes',[360-72 360-66],'latitudes',[12,18]);    hold on
    for i = 1:length(k)
        m_plot(k(i).edgelon,k(i).edgelat,'-','linewi',1,'color',[.4 .3 .7]);
    end
    m_coast('patch',[0.5 0.5 0.5],'edgecolor',[0.2 0.2 0.2]);
    m_grid('xtick',5,'box','fancy','tickdir','out','LineWidth',0.5);
%     print(gcf,'-djpeg','-r450',['E:\Figure\','Contour_surf.jpg']);
    %% Submesoscale eddies
    eddy_contour = get_ow_eddy(lon_m, lat_m, W_noiseless, zeta, OWValue, ShapeError); % 提取涡旋轮廓
    figure%('visible','off')
    set(gcf,'outerposition',get(0,'screensize'),'color','w');
    set(gca,'fontsize',26,'Fontname','Time new roman','fontweight','bold');hold on; 
    m_proj('miller','longitudes',[360-72 360-66],'latitudes',[12,18]);    hold on
    for i = 1:length(eddy_contour)
        m_plot(eddy_contour(i).edgelon,eddy_contour(i).edgelat,'-','linewi',1,'color',[.5 .5 .5]);
    end
    m_coast('patch',[0.5 0.5 0.5],'edgecolor',[0.2 0.2 0.2]);
    m_grid('xtick',5,'box','fancy','tickdir','out','LineWidth',0.5);
%     print(gcf,'-djpeg','-r450',['E:\Figure\','Eddy_Contour_surf.jpg']);