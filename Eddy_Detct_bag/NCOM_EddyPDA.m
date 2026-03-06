function PDA = NCOM_EddyPDA(EddyData, Eddymodel, Count_num)
% 计算NCOM 3D涡旋的（垂向）位势密度异常PDA（0~300m深度变化）
%   输入：3D涡旋的涡旋轮廓EddyData，网格Eddymodel，追踪Count信息
%   输出：各层级的涡内均值PDA
    field = fieldnames(EddyData); 
    deplayer = [1,3,6,9,11,13,15,16,18,20:25];
    LL = Count_num.layer(:,1); % 层
    index = Count_num.layer(:,2); % 索引
    PDA = []; 
    for dd = 1:15 %size(Count(num).layer,1)
        % 读取涡旋（主要提取密度）
        if ismember(dd,LL) % ----------若某一层“存在”追踪的涡旋----------
            eddy_para = EddyData.(field{dd})(index(dd==LL));
            eddylon = eddy_para.lon;
            eddylat = eddy_para.lat;
            eddy_density = eddy_para.potential_density_hp;
            edgelon = eddy_para.edgelon;
            edgelat = eddy_para.edgelat;
        elseif dd<LL(1) % --若某一层“不存在”追踪的涡旋，且层数在涡旋上---
            eddylon = Eddymodel.lon;
            eddylat = Eddymodel.lat;
            eddy_density = Eddymodel.potential_density_hp(:,:,deplayer(dd));
            edgelon = EddyData.(field{LL(1)})(index(1)).edgelon;
            edgelat = EddyData.(field{LL(1)})(index(1)).edgelat;
        else% --------------若某一层“不存在”追踪的涡旋，且层数在涡旋下---
            eddylon = Eddymodel.lon;
            eddylat = Eddymodel.lat;
            eddy_density = Eddymodel.potential_density_hp(:,:,deplayer(dd));
            edgelon = EddyData.(field{LL(end)})(index(end)).edgelon;
            edgelat = EddyData.(field{LL(end)})(index(end)).edgelat;
        end
        [xlen,ylen] = size(eddy_density);
        [lat_f,lon_f] = meshgrid(eddylat,eddylon);
        % 层级参数
        res =0.05; % 提升分辨率（原始分辨率太低，求出的涡内均值不准确）
        lon = interp1(1:xlen,eddylon,1:res:xlen);
        lat = interp1(1:ylen,eddylat,1:res:ylen);
        [lat_re,lon_re] = meshgrid(lat,lon);
        density = interp2(lat_f,lon_f,eddy_density,lat_re,lon_re,'linear');
        clear eddylon eddylat eddy_density xlen ylen lon lat lat_f lon_f
        % ---涡旋内---
        in=inpolygon(lon_re, lat_re, edgelon, edgelat); %内部点
        density(~in)=nan;
        PDA(dd,1) = mean(density,'all','omitnan');        
        clear in density lat_re lon_re edgelon edgelat eddy_para
    end
end

