function k = get_ow_eddy(lon, lat, ow, vort, varargin)
% 筛选OW参数内部值和形状误差符合要求的轮廓，为涡旋轮廓
% lon, lat：输入要求和ow吻合的二维经纬度网格
% ow：提取涡旋轮廓
% vort：判断涡旋极性
    p = inputParser;
    p.addRequired('lon',@(x)validateattributes(x,{'numeric'},...
        {'nonempty'},'eddyDetection_output','lon',1));
    p.addRequired('lat',@(x)validateattributes(x,{'numeric'},...
        {'nonempty'},'eddyDetection_output','lat',2));
    p.addRequired('ow',@(x)validateattributes(x,{'numeric'},...
        {'3d'},'eddyDetection_output','watervar',3));
    p.addRequired('vort',@(x)validateattributes(x,{'numeric'},...
        {'3d'},'eddyDetection_output','watervar',4));
% -----------------------------varargin------------------------------------
    defaultValue=-0.2; % 默认直径，unit: Km
    p.addOptional('OWValue',defaultValue,...
        @(x)validateattributes(x,{'numeric'},...
        {'nonzero'},'eddyDetection_output','radiusThreshold',5));
    defaultShapeError=0.35; % 默认形状误差
    p.addOptional('ShapeError',defaultShapeError,...
        @(x)validateattributes(x,{'numeric'},...
        {'nonzero','positive'},'eddyDetection_output','ShapeError',6));
    p.parse(lon,lat,ow,vort,varargin{:});
    
    Value=p.Results.OWValue;%OW阈值
    Error=p.Results.ShapeError; % 形状误差
%--------------------------------------------------------------------------
    figure('visible','off');
    [c,~]=contour(lon, lat, ow, [Value 9999]);
    s = getcontourlines(c); % getcontourlines 解析contour得到的"闭合"轮廓
    ss=0;
    k = struct();
    warning off
    for i = 1:length(s) 
    % 保留内部数值<-0.2的轮廓线(Williams et al. 2011)
        in=inpolygon(lon, lat, s(i).x, s(i).y); %内部点
        if sum(in,'all')>=5 %内部至少存在5个像素点
            if mean(ow(in),'omitnan') <= Value  %内部值<-0.2
                % 保留近圆形的轮廓（形状误差<35%）(Kurian et al. 2011)
                x = s(i).x;   %经度
                y = s(i).y;   %纬度
                center = [mean(x),mean(y)]; % 拟合圆心
                [arclen,az]=distance(y,x,center(2),center(1));
                dis = arclen./180.*pi.*6371; % 单位：km
                [xx,yy]=pol2cart(az*pi/180,dis);
                R= mean(dis); % 拟合圆半径
                A = pi.*R.^2; % 拟合圆面积
                A2 = polyarea(xx,yy); % 多边形面积
                inv = pi/25; % 角度间隔
                [dx,dy]=pol2cart([0:inv:2*pi-inv,0],repmat(R,1,51)); % 50个点的近似圆形（2pi/50）
                polyout = intersect(polyshape(xx,yy),polyshape(dx,dy) );
                shape_error = (A+A2-2*polyout.area)./A;
                if shape_error < Error % 形状误差<35%
                    ss = ss+1;
                    if mean(vort(in),'omitnan') > 0 % 逆时针为+，顺时针为-
                        k(ss).type='C'; % 气旋
                    else
                        k(ss).type='A'; % 反气旋
                    end
                    k(ss).center = center;
                    k(ss).edgelon = s(i).x;
                    k(ss).edgelat = s(i).y;
                    k(ss).radius = R;
                    k(ss).shape_error = shape_error;
                end
            end
        end
    end
    close all
end
