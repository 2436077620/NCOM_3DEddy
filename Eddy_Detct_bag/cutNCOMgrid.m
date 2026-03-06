function Eddymodel = cutNCOMgrid(Model,Cutlim)
% 裁剪读取的Model文件（读取的子函数：ReadNCOMmodel.m）
    % 输入：Model、Cutlim
    % ----------------模式数据与Eddy Cutlim匹配，裁剪----------------------
    model_field = fieldnames(Model);
    [~,lon_i] = min(abs(Model.lon-Cutlim(1,:)));
    [~,lat_i] = min(abs(Model.lat-Cutlim(2,:)));
    if lon_i(1) > 1
        lon_i(1)=lon_i(1)-1; 
    end
    if lon_i(2) < length(Model.lon)
        lon_i(2)=lon_i(2)+1;
    end
    if lat_i(1) > 1
     lat_i(1)=lat_i(1)-1; 
    end
    if lat_i(2) < length(Model.lat)
        lat_i(2)=lat_i(2)+1;
    end
    Eddymodel.lon = Model.lon(lon_i(1):lon_i(2)); % 经度
    Eddymodel.lat = Model.lat(lat_i(1):lat_i(2)); % 纬度
    Eddymodel.depth = Model.depth; % 深度
    [lat_xy,~] = meshgrid(Eddymodel.lat,Eddymodel.lon);
    % 水下密度异常3D网格数据
    for ff=5:length(model_field) % 提取Output网格参数
        Eddymodel.(model_field{ff}) = Model.(model_field{ff})(lon_i(1):lon_i(2),lat_i(1):lat_i(2),:);
    end   
    % ----相对涡度----
    [~,dudy] = gradient(permute(Eddymodel.water_u_hp,[2,1,3]),(1/30)/180*pi*6371*10^3); % u
    [dvdx,~] = gradient(permute(Eddymodel.water_v_hp,[2,1,3]),(1/30)/180*pi*6371*10^3); % v
    vorticity = permute(dvdx-dudy,[2,1,3]); % 逆时针为+，顺时针为-
    Eddymodel.zeta = vorticity./gsw_f(lat_xy);% 利用f归一化    
end

