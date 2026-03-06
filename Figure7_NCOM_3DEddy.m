clc;clear;close all;
% 选择样本，绘制垂向追踪后的3D Example （提取出的亚中尺度涡三维结构）
% 包括 1）垂向追踪涡旋轮廓 2）提取的涡旋内部参数
% ============统计track个数================
load('E:\DATA\NCOMEddy\Eddy_coamps_ncom_amseas_u_1_2023072600_00210000.mat');
field = fieldnames(EddyData);
clear track 
for dd = 1:length(field) % 层级
    if dd==1
        track = [EddyData.(field{dd}).track]';
    else
        track = [track;[EddyData.(field{dd}).track]'];
    end
end
count = tabulate(track);
track_num = 304;%track_num(3);
%% 画图
posit = [700, 0.63, 0.82]; % gcf宽度，gca右边界，colorbar左边界
var_all = {'vorticity','uv', 'TA', 'SA', 'PDA'};
for vv = 1:5
    var_name = var_all{vv};
    if strcmp(var_name, 'vorticity')
        var_str = '\zeta/f';
        var_range= [-0.5,0.5];
    elseif strcmp(var_name, 'TA')
        var_str = 'TA (\circC)';
        var_range= [-0.2,0.21];
    elseif strcmp(var_name, 'SA')
        var_str = 'SA (\times10^-^2 psu)';
        var_range= [-2,2];
    elseif strcmp(var_name, 'PDA')
        var_str = 'PDA (\times10^-^2 kg/m^3)';
        var_range= [-6,6];
    elseif strcmp(var_name, 'w')
        var_str = 'w (10^-^3 m/s)';
        var_range= [-1.5,1.5];
    elseif strcmp(var_name, 'uv')
        var_str = 'Velocity (m/s)';
        var_range= [0,0.12];
    end
    % ================================
    figure;
    set(gcf,'position',[200,50,posit(1),950],'color','w');
    set(gca,'position',[0.17,0.08,posit(2),0.9],'FontName','Times New Roman','FontSize',22,'FontWeight','bold');hold on;
    lonlim = []; latlim = []; eddycen = [];
    %--------------绘制参数-----------------
    for dd=1:length(field)
        track1 = [EddyData.(field{dd}).track]';
        index = find(track1 == track_num);
        if isempty(index)
            continue;
        end
        dep = -str2num(field{dd}(7:end-1)); %深度
        eddy_para = EddyData.(field{dd})(index); % 被追踪的涡旋
        [X,Y] = LL3XY(eddy_para.center(1),eddy_para.center(2));
        eddycen(dd,:) = [X, Y, dep]; 
        if ~ismember(dd,[1,4,6,8:length(field)]) % 仅绘制特定层的图像
            continue;
        end
        lonlim = [min([eddy_para.edgelon, lonlim]),max([eddy_para.edgelon, lonlim])]; % 范围
        latlim = [min([eddy_para.edgelat, latlim]),max([eddy_para.edgelat, latlim])]; % 范围
    %     deplim = [dep,max(dep,deplim)];
        [lat_f,lon_f] = meshgrid(eddy_para.lat,eddy_para.lon);
        if strcmp(var_name, 'vorticity')
            var0 = eddy_para.vorticity./gsw_f(lat_f); % 绘制的参数（相对涡度需要基于f归一化）
        elseif strcmp(var_name, 'TA')
            var0 = eddy_para.temp_hp;
        elseif strcmp(var_name, 'SA')
            var0 = eddy_para.salinity_hp.*100;
        elseif strcmp(var_name, 'PDA')
            var0 = eddy_para.potential_density_hp.*100;
        elseif strcmp(var_name, 'uv')
            ua = eddy_para.water_u_hp;
            va = eddy_para.water_v_hp;
            w = eddy_para.water_w;
            var0 = hypot(ua, va);
        end
        [xlen,ylen] = size(var0);
        % ----插值----
        % 层级参数
        res =0.02; % 提升分辨率
        lon = interp1(1:xlen,eddy_para.lon,1:res:xlen);
        lat = interp1(1:ylen,eddy_para.lat,1:res:ylen);
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
        var(~in)=nan;
        % ----画图----
        [x1,y1] = LL3XY(lon_re,lat_re);
        [x2,y2] = LL3XY(lon_f, lat_f);
        surf(x1-eddycen(1,1),y1-eddycen(1,2),dep.*ones(size(lat_re)),var); shading interp% 在z=h的位置绘制图层
        [e1,e2] = LL3XY(edgelon,edgelat);
        if strcmp(var_name, 'uv')
            in0 = inpolygon(lon_f, lat_f, edgelon, edgelat);
            ua(~in0)=nan; va(~in0)=nan; w(~in0)=nan;
            quiver3(x2-eddycen(1,1),y2-eddycen(1,2), dep.*ones(size(lat_f)),ua,va,w,'color',[.6 .6 .6],'linewidth',1);
        else
           plot3(e1-eddycen(1,1), e2-eddycen(1,2), repmat(dep,1,(edgelen-1)/res2+1),'-k','linewi',0.6,'color',[.4 .4 .4]);
        end
        clear dep lat_re lon_re var edgelen edgelon edgelat in lat_f lon_f in0
    end
    if strcmp(var_name, 'uv')
        colormap(flip(othercolor('RdYlBu3'))); %Spectral10')));% 
    else
        colormap(m_colmap('diverging'));
    end
    h1=colorbar('location','east','position',[posit(3) 0.14 0.025 0.68]);
    h1.Label.String = var_str;
    caxis(var_range);
    % 涡旋中心的垂向连线
    eddycen0 = [eddycen(:,1)-eddycen(1,1), eddycen(:,2)-eddycen(1,2)];
    plot3(eddycen0(:,1), eddycen0(:,2), eddycen(:,3), '--' ,'linewi',1.5, 'color', [102 109 101]./255);
    plot3(eddycen0(:,1), eddycen0(:,2), eddycen(:,3)+0.2, '.','MarkerSize',20,'color', [102 109 101]./255)
    view(-28,11);
    set(gca,'xlim',[-12.6,13],'ylim',[-12.5,17.2],'zlim',[eddycen(end,3)-0.1,0.2],'color','w')
    pos = axis;
    xlabel('dX (km)','rotation',12, 'position',[pos(1)+0.35*(pos(2)-pos(1)), pos(3)-0.32*(pos(4)-pos(3)), pos(5)]);
    ylabel('dY (km)', 'rotation',-45, 'position',[pos(1)-0.32*(pos(2)-pos(1)), pos(3)+0.12*(pos(4)-pos(3)), pos(5)]);
    zlabel('Depth (m)');
    box on
    grid on
    % print(gcf,'-djpeg','-r450',['E:\Figure\Example_',var_name,'_304.jpg']);
end