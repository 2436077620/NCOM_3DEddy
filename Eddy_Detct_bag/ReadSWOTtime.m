function swottime = ReadSWOTtime(name)
    % SWOT从数据名称中提取时间（比直接读取更快）
    year =  [str2num(name(:,31:34)), str2num(name(:,47:50))];
    month = [str2num(name(:,35:36)), str2num(name(:,51:52))];
    day =   [str2num(name(:,37:38)), str2num(name(:,53:54))];
    hour =  [str2num(name(:,40:41)), str2num(name(:,56:57))];
    min =   [str2num(name(:,42:43)), str2num(name(:,58:59))];
    sec =   [str2num(name(:,44:45)), str2num(name(:,60:61))];
    swottime = (datenum(year(:,1),month(:,1),day(:,1),hour(:,1),min(:,1),sec(:,1))+...
        datenum(year(:,2),month(:,2),day(:,2),hour(:,2),min(:,2),sec(:,2)))./2; %起始和终止时间的均值
end

