function sst = ReadSSTfile(pathname,cutlon,cutlat)
% 基于目标时间（swottime），读取对应时间范围的SSS文件夹（文件）
    sst = struct(); 
    
    % 读所有文件及时间范围
    sstfile  = dir(pathname);% 所有SST文件
    file = [sstfile.folder,'\',sstfile.name];
    lon = double(ncread(file,'lon'));
    lat = double(ncread(file,'lat'));
    
    %读取当天的SST & 裁剪
    [~,index_lon] = min(abs(lon-cutlon));
    [~,index_lat] = min(abs(lat-cutlat));
    start = [index_lon(1)-1, index_lat(1)-1,1];
    len = [index_lon(2)-index_lon(1)+3, index_lat(2)-index_lat(1)+3,1];
    sst.temp = ncread(file,'analysed_sst',start,len);
    sst.tempa = ncread(file,'sst_anomaly',start,len);
    
    [sst.lat,sst.lon] = meshgrid(lat(index_lat(1)-1:index_lat(2)+1), lon(index_lon(1)-1:index_lon(2)+1));
end