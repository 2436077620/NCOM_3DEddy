function swotGrid = cutswotEddy(swot,swotEddy)
    % 裁剪检测出的SWOT Eddy的海面高网格（用于主函数SWOTEddyDetct.m）
    % 输入：已经读取的swot原始文件swot
    % 输入：检测出的涡旋：swotEddy
    % 输出：裁剪并插值到规则网格中的swot数据swotGrid
    swotGrid = struct();
    Eddycenlon = swotEddy.center(1);
    Eddycenlat = swotEddy.center(2);
    LonDis = max(swotEddy.edge(1,:))-min(swotEddy.edge(1,:));% 经度差 与所在纬度有关
    LatDis = max(swotEddy.edge(2,:))-min(swotEddy.edge(2,:));% 纬度差 直接转为弧度
    % 计算Eddy 2倍空间范围， 建立 "裁剪框多边形"
    cutlat = [Eddycenlat-LatDis, Eddycenlat+LatDis, Eddycenlat+LatDis, Eddycenlat-LatDis, Eddycenlat-LatDis];
    cutlon = [Eddycenlon-LonDis, Eddycenlon-LonDis, Eddycenlon+LonDis, Eddycenlon+LonDis, Eddycenlon-LonDis];
    % 裁剪SWOT变量
    [swot_in,~] = inpolygon(swot.lon,swot.lat,cutlon,cutlat);
    lon = swot.lon(swot_in);
    lat = swot.lat(swot_in);
    % 裁剪数据插值到2D规则网格，空间分辨率0.015°
%     swotGrid.time = swot.time;
    swotGrid.lon = [round(min(lon),2):0.015:round(max(lon),2)]';
    swotGrid.lat = [round(min(lat),2):0.015:round(max(lat),2)]';
    [lat_xy,lon_xy]=meshgrid(swotGrid.lat, swotGrid.lon);
    fields = fieldnames(swot);
    for ff = 4:length(fields)
        swotGrid.(fields{ff}) = griddata(lon, lat, swot.(fields{ff})(swot_in), lon_xy, lat_xy);
    end
end

