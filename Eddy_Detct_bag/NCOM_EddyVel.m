function Velocity = NCOM_EddyVel(EddyData, Eddymodel_num, Count_num)
% 计算NCOM 3D涡旋的（垂向）水平流场及相关参数（0~300m深度变化）
% 背景流、KE、KE强度
% 亚中尺度涡流场、EKE、EI、涡度
%   输入：3D涡旋的涡旋轮廓EddyData，网格Eddymodel，追踪Count信息
%   就涡度参数而言，在Eddymodel被存为 zeta=vorticity./gsw_f(lat);而EddyData没有./f
%   输出：各层级的涡内均值PDA
    % 提取垂向参数
    field = fieldnames(EddyData); 
    deplayer = [1,3,6,9,11,13,15,16,18,20:25];
    LL = Count_num.layer(:,1); % 层
    index = Count_num.layer(:,2); % 索引
    Velocity = struct();
    for dd = 1:15 %size(Count(num).layer,1)
        % 读取涡旋（主要提取密度）
        if ismember(dd,LL) % ----------若某一层“存在”追踪的涡旋----------
            eddy_para = EddyData.(field{dd})(index(dd==LL));
            eddylon = eddy_para.lon;
            eddylat = eddy_para.lat;
            [lat_f,lon_f] = meshgrid(eddylat,eddylon);
            % 原始流场
            u0 = eddy_para.water_u;
            v0 = eddy_para.water_v;            
            % 亚中尺度流参数
            ua0=eddy_para.water_u_hp; 
            va0=eddy_para.water_v_hp; 
            EKE0=(ua0.^2+va0.^2)./2; 
            zeta0 = eddy_para.vorticity./gsw_f(lat_f);
            % 背景流参数（原始流场 - 高通流场）
            ub0=u0-ua0; 
            vb0=v0-va0; 
            KE0 = (ub0.^2+vb0.^2)./2;
            edgelon = eddy_para.edgelon;
            edgelat = eddy_para.edgelat;
            R = eddy_para.radius;            
        elseif dd<LL(1) % --若某一层“不存在”追踪的涡旋，且层数在涡旋上---
            eddylon = Eddymodel_num.lon;
            eddylat = Eddymodel_num.lat;
            [lat_f,lon_f] = meshgrid(eddylat,eddylon);
            % 原始流场
            u0 = Eddymodel_num.water_u(:,:,deplayer(dd));
            v0 = Eddymodel_num.water_v(:,:,deplayer(dd));            
            % 亚中尺度流参数
            ua0=Eddymodel_num.water_u_hp(:,:,deplayer(dd)); 
            va0=Eddymodel_num.water_v_hp(:,:,deplayer(dd)); 
            EKE0=(ua0.^2+va0.^2)./2; 
            zeta0 = Eddymodel_num.zeta(:,:,deplayer(dd));
            % 背景流参数
            ub0=u0-ua0; 
            vb0=v0-va0; 
            KE0 = (ub0.^2+vb0.^2)./2;
            edgelon = EddyData.(field{LL(1)})(index(1)).edgelon;
            edgelat = EddyData.(field{LL(1)})(index(1)).edgelat;
            R = EddyData.(field{LL(1)})(index(1)).radius;
        else% --------------若某一层“不存在”追踪的涡旋，且层数在涡旋下---
            eddylon = Eddymodel_num.lon;
            eddylat = Eddymodel_num.lat;
            [lat_f,lon_f] = meshgrid(eddylat,eddylon);
            % 原始流场
            u0 = Eddymodel_num.water_u(:,:,deplayer(dd));
            v0 = Eddymodel_num.water_v(:,:,deplayer(dd));            
            % 亚中尺度流参数
            ua0=Eddymodel_num.water_u_hp(:,:,deplayer(dd)); 
            va0=Eddymodel_num.water_v_hp(:,:,deplayer(dd)); 
            EKE0=(ua0.^2+va0.^2)./2; 
            zeta0 = Eddymodel_num.zeta(:,:,deplayer(dd));
            % 背景流参数
            ub0=u0-ua0; 
            vb0=v0-va0; 
            KE0 = (ub0.^2+vb0.^2)./2;
            edgelon = EddyData.(field{LL(end)})(index(end)).edgelon;
            edgelat = EddyData.(field{LL(end)})(index(end)).edgelat;
            R = EddyData.(field{LL(end)})(index(end)).radius;
        end
        [xlen,ylen] = size(u0);
        % 层级参数
        res =0.05; % 提升分辨率（原始分辨率太低，求出的涡内均值不准确）
        lon = interp1(1:xlen,eddylon,1:res:xlen);
        lat = interp1(1:ylen,eddylat,1:res:ylen);
        [lat_re,lon_re] = meshgrid(lat,lon);
        % 亚中尺度流参数
        ua = interp2(lat_f,lon_f,ua0,lat_re,lon_re,'linear');
        va = interp2(lat_f,lon_f,va0,lat_re,lon_re,'linear');
        EKE = interp2(lat_f,lon_f,EKE0,lat_re,lon_re,'linear');
        zeta = interp2(lat_f,lon_f,zeta0,lat_re,lon_re,'linear');
        % 背景流参数
        ub = interp2(lat_f,lon_f,ub0,lat_re,lon_re,'linear');
        vb = interp2(lat_f,lon_f,vb0,lat_re,lon_re,'linear');
        KE = interp2(lat_f,lon_f,KE0,lat_re,lon_re,'linear');
        % ---涡旋内---
        in = inpolygon(lon_re, lat_re, edgelon, edgelat); %内部点
        EKE(~in) = nan;  zeta(~in) = nan; 
        ub(~in) = nan; vb(~in) = nan; KE(~in) = nan;
        % ---涡边缘---
        index_e=false(size(ua));
        for ii=1:length(edgelon)-1
            [~,index_x] = min(abs(lon-edgelon(ii)),[],'omitnan');
            [~,index_y] = min(abs(lat-edgelat(ii)),[],'omitnan');
            index_e(index_x,index_y)=true;
        end
        edge_velocity = mean(hypot(ua(index_e),va(index_e)),'all','omitnan');
        clear eddylon eddylat eddy_density xlen ylen lon lat lat_f lon_f
        % 亚中尺度参数        
        Velocity.velocity_edge(dd,1) = edge_velocity; % 平均边缘流速
        Velocity.EKE(dd,1) = sum(EKE,'all','omitnan'); % 总的涡动能
        Velocity.EI(dd,1) = Velocity.EKE(dd,1)./(pi*R^2);% 涡旋强度（涡旋单位面积上的EKE）
        Velocity.zeta(dd,1) = mean(zeta,'all','omitnan'); % 平均涡度
        % 背景流（大尺度）参数
        Velocity.ub(dd,1) = mean(ub,'all','omitnan');
        Velocity.vb(dd,1) = mean(vb,'all','omitnan');
        Velocity.KE(dd,1) = sum(KE,'all','omitnan');
        clear in density lat_re lon_re edgelon edgelat eddy_para
    end
end

