clc;clear;close all;
% 研究区的涡度、高通滤波后的涡度
region0  = {'alaska','amseas'}; % 涡旋所在海域  alaska, amseas
eddypath = 'E:\DATA\NCOMGrid\';
res = 3.3;
hpkm = 50;  % 对模式数据高通滤波的尺寸
for i = 1:2
    region = region0{i};
    model_path= ['E:\DATA\NCOM_model\',region,'\'];     % 储存模式数据的路径
    % ---------------------------NCOM模式数据的范围----------------------------
    if strcmp(region,'alaska')
        eddy_file = '_008_498_';
        range = [195, 238, 36.45, 61]; % 0~360
        box = [360-149 360-143 54,58.5];
    elseif strcmp(region,'amseas')
        eddy_file = '_001_173_'; 
        range = [261.9, 305.1, 5, 31.1];  
        box = [360-72 360-66 12,18];
    end
    figrange = [range(1:2);range(3:4)];
    eddydata = dir([eddypath,'*',eddy_file,'*.mat']);
    load([eddypath,eddydata.name]);
    Model = ReadNCOMmodel([model_path, ModelData.NCOM_Name],res,hpkm);% 读model文件
    %% NCOM Vorticity & Okubo-Weiss 参数计算
    [lat_m, lon_m] = meshgrid(Model.lat, Model.lon);
    surf_el_hp = filt2(Model.surf_el,res,hpkm,'hp');
    layer = 1;
    % --------------------原始涡度----------------------------
    u = Model.water_u(:,:,layer); 
    v = Model.water_v(:,:,layer);
    [~,dudy] = gradient(permute(u,[2,1,3]),(1/30)/180*pi*6371*10^3); % u
    [dvdx,~] = gradient(permute(v,[2,1,3]),(1/30)/180*pi*6371*10^3); % v
    vorticity = permute(dvdx-dudy,[2,1,3]); % 逆时针为+，顺时针为-
    zeta = vorticity./gsw_f(lat_m);% 利用f归一化   
    clear dudx dudy dvdx dvdy
    % --------------------filtered涡度----------------------------
    u_hp = filt2(Model.water_u(:,:,layer),res,hpkm,'hp');
    v_hp = filt2(Model.water_v(:,:,layer),res,hpkm,'hp');
    [~,dudy] = gradient(permute(u_hp,[2,1,3]),(1/30)/180*pi*6371*10^3); % u
    [dvdx,~] = gradient(permute(v_hp,[2,1,3]),(1/30)/180*pi*6371*10^3); % v
    vorticity = permute(dvdx-dudy,[2,1,3]); % 逆时针为+，顺时针为-
    zeta_hp = vorticity./gsw_f(lat_m);% 利用f归一化   
    clear dudx dudy dvdx dvdy
    %%  --------------(原始)涡度------------------------
    figure%('visible','off')%-----------------------------------------
    set(gcf,'outerposition',get(0,'screensize'),'color','w');
    set(gca,'fontsize',26,'Fontname','Time new roman','fontweight','bold');hold on; 
    m_proj('miller','longitudes',figrange(1,:),'latitudes',figrange(2,:));
    hold on
    m_pcolor(lon_m, lat_m, zeta);shading interp;
    colormap(flip(othercolor('RdBu11')));
    c = colorbar('fontsize',26,'location','southout');
    c.Label.String = '\zeta/f';
    if strcmp(region,'alaska')
        caxis([-0.301 0.3]);
    elseif strcmp(region,'amseas')
        caxis([-0.801 0.8]);
    end
    m_gshhs_i('patch',[0.5 0.5 0.5],'edgecolor',[0.2 0.2 0.2]);
    m_plot([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)],'-k','linewi',3);
    m_grid('xtick',5,'box','fancy','tickdir','out','LineWidth',0.5);
%     print(gcf,'-djpeg','-r450',['E:\Figure\','Vort_', num2str(Model.depth(layer)),'m_raw_all.jpg']); 
    %% --------------(滤波后的)涡度------------------------
    figure%('visible','off')%-----------------------------------------
    set(gcf,'outerposition',get(0,'screensize'),'color','w');
    set(gca,'fontsize',26,'Fontname','Time new roman','fontweight','bold');hold on; 
    m_proj('miller','longitudes',figrange(1,:),'latitudes',figrange(2,:));
    m_pcolor(lon_m, lat_m, zeta_hp);shading interp;
    colormap(flip(othercolor('RdBu11')));
    c = colorbar('fontsize',26,'location','southout');
    c.Label.String = '\zeta/f';
    if strcmp(region,'alaska')
        caxis([-0.151 0.15])
    elseif strcmp(region,'amseas')
        caxis([-0.401 0.4]);
    end
    m_gshhs_i('patch',[0.5 0.5 0.5],'edgecolor',[0.2 0.2 0.2]);
    m_plot([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)],'-k','linewi',3);
    m_grid('xtick',5,'box','fancy','tickdir','out','LineWidth',0.5);
%     print(gcf,'-djpeg','-r450',['E:\Figure\','Vort_', num2str(Model.depth(layer)),'m_50km_all.jpg']); 
    %%  --------------子区域(原始)涡度------------------------
    figure%('visible','off')
    set(gcf,'outerposition',get(0,'screensize'),'color','w');
    set(gca,'fontsize',26,'Fontname','Time new roman','fontweight','bold');hold on; 
    m_proj('miller','longitudes',[box(1) box(2)],'latitudes',[box(3),box(4)]);
    hold on
    m_pcolor(lon_m, lat_m, zeta);shading interp;
    colormap(flip(othercolor('RdBu11')));
    if strcmp(region,'alaska')
        caxis([-0.301 0.3]);
    elseif strcmp(region,'amseas')
        caxis([-0.801 0.8]);
    end
    m_gshhs_i('patch',[0.5 0.5 0.5],'edgecolor',[0.2 0.2 0.2]);
    m_plot([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)],'-k','linewi',3);
    m_grid('xtick',5,'xticklabel',[],'yticklabel',[],'box','on','tickdir','in','LineWidth',0.5);
%     print(gcf,'-djpeg','-r450',['E:\Figure\','Vort_', num2str(Model.depth(layer)),'m_raw_1.jpg']); 
    %% --------------子区域(滤波后的)涡度------------------------
    figure%('visible','off')
    set(gcf,'outerposition',get(0,'screensize'),'color','w');
    set(gca,'fontsize',26,'Fontname','Time new roman','fontweight','bold');hold on; 
    m_proj('miller','longitudes',[box(1) box(2)],'latitudes',[box(3),box(4)]);    hold on
    m_pcolor(lon_m, lat_m, zeta_hp);shading interp;
    colormap(flip(othercolor('RdBu11')));
    if strcmp(region,'alaska')
        caxis([-0.151 0.15])
    elseif strcmp(region,'amseas')
        caxis([-0.401 0.4]);
    end
    m_gshhs_i('patch',[0.5 0.5 0.5],'edgecolor',[0.2 0.2 0.2]);
    m_plot([box(1) box(1) box(2) box(2) box(1)], [box(3) box(4) box(4) box(3) box(3)],'-k','linewi',3);
    m_grid('xtick',5,'xticklabel',[],'yticklabel',[],'box','on','tickdir','in','LineWidth',0.5);
%     print(gcf,'-djpeg','-r450',['E:\Figure\','Vort_', num2str(Model.depth(layer)),'m_50km_1.jpg']); 
end