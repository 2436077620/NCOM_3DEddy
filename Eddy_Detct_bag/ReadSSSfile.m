function sss = ReadSSSfile(pathname,swottime)
% 基于目标时间（swottime），读取对应时间范围的SSS文件夹（文件）
    sss = struct();
    
    % 读所有文件及时间范围
    sssfile0  = dir(pathname);% 所有SSS文件
    sss_timerange = zeros(length(sssfile0),2);
    for i = 1:length(sssfile0)
        num = regexp(sssfile0(i).name,'\d+','match');
        sss_timerange(i,:) =  [datenum(num{2},'yyyymmdd'), datenum(num{3},'yyyymmdd')];
    end
    index = find(swottime>=sss_timerange(:,1) & swottime<sss_timerange(:,2)+1);
    file = [sssfile0(index).folder, '\', sssfile0(index).name];
    
    time = ncread(file,'time')./86400 + datenum(1970,1,1);
    lon = double(ncread(file,'longitude'));
    lat = double(ncread(file,'latitude'));
    [sss.lat,sss.lon] = meshgrid(lat,lon);
    
    %读取当天的SSS
    [~,index_time] = min(abs(time-swottime));
    start = [1, 1, 1, index_time];
    len = [length(lon), length(lat), 1, 1];
    sss.salt = ncread(file,'sos',start,len);
end

