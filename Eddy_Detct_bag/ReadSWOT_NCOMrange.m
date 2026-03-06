function swot = ReadSWOT_NCOMrange(swotfile,range)
    % 读取经过特定空间范围的SWOT文件（用于主函数SWOTEddyDetct.m）
    % 读取数据版本v1.0.2、v2.0.1
    version = swotfile(end-7:end-3);
    % 输入swot路径文件，range为指定范围
    swotlon=double(ncread(swotfile,'longitude')); %0~360
    swotlat=double(ncread(swotfile,'latitude')); 
    % 判断SWOT是否经过指定范围
    [x1,y1] = find(swotlon>=range(1) & swotlon<=range(2) & swotlat>=range(3) & swotlat<=range(4));
    swot = struct();%若不在，则输出空结构体
    if ~isempty(x1) && ~isempty(y1) % 若在，则裁剪指定范围
            x=min(x1); y=min(y1); dx=max(x1)-x+1; dy=max(y1)-y+1;
            swot.time = mean(double(ncread(swotfile,'time')))/86400 + datenum(2000,1,1);
            swot.lon=double(ncread(swotfile,'longitude',[x,y],[dx,dy]));  %0~360
            swot.lat=double(ncread(swotfile,'latitude',[x,y],[dx,dy])); 
            swot.mdt = double(ncread(swotfile,'mdt',[x,y],[dx,dy]));      % MDT
            swot.sigma0 = double(ncread(swotfile,'sigma0',[x,y],[dx,dy]));
        if strcmp(version,'1.0.2') || strcmp(version,'_v1.0')
            swot.ssha=double(ncread(swotfile,'ssha_noiseless',[x,y],[dx,dy])); %SSHA without noise
            swot.ug = double(ncread(swotfile,'ugos',[x,y],[dx,dy]));
            swot.vg = double(ncread(swotfile,'vgos',[x,y],[dx,dy]));
            swot.ua = double(ncread(swotfile,'ugosa',[x,y],[dx,dy]));
            swot.va = double(ncread(swotfile,'vgosa',[x,y],[dx,dy]));
        elseif strcmp(version,'2.0.1') 
            swot.ssha=double(ncread(swotfile,'ssha_filtered',[x,y],[dx,dy])); %SSHA without noise
            swot.ug = double(ncread(swotfile,'ugos_filtered',[x,y],[dx,dy]));
            swot.vg = double(ncread(swotfile,'vgos_filtered',[x,y],[dx,dy]));
            swot.ua = double(ncread(swotfile,'ugosa_filtered',[x,y],[dx,dy]));
            swot.va = double(ncread(swotfile,'vgosa_filtered',[x,y],[dx,dy]));
        end
    end
end

