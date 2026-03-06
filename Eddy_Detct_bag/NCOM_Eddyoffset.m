function [direct,angle] = NCOM_Eddyoffset(Count_num,Depth)
% 计算NCOM 3D涡旋的（垂向）偏移：偏移方位（相较于正北）&偏移角度（拟合线与垂线夹角）
% 内含自定义函数 uv_angle
    center = Count_num.center; % 中心位置
    [X,Y]=LL3XY(center(:,1),center(:,2)); %（经纬度->直角坐标）输出单位km
    X = X*1000; Y = Y*1000; % 单位变成：m
    % ----------插值到均匀深度-----------
    depth = Depth(Count_num.layer(:,1)); % 深度
    z = (depth(1):-2:depth(end))';          % 插值到2m间隔
    x = interp1(depth, X, z, 'linear');
    y = interp1(depth, Y, z, 'linear');
    %         z = depth;x=X;y=Y;
    % 计算平均值（拟合的直线必过所有坐标的算数平均值）
    x = x-mean(x);  y = y-mean(y);% 拟合点坐标（均值放到原点）
    % ========奇异值分解计算方向向量=========
    % 协方差矩阵奇异变换：所得直线的方向实际上与最大奇异值对应的奇异向量相同
    centeredLine=bsxfun(@minus,[x,y,z],[0 0 mean(z)]); % 经过原点
    [~,~,V]=svd(centeredLine);
    direction=V(:,1);%方向向量
    % 三维线性回归
    zz = -300; % 300m深度
    t = (zz-mean(z))./direction(3); 
    xx = direction(1)*t;
    yy = direction(2)*t;
    direct = uv_angle(xx,yy); % uv_angle为自定义函数
    angle = atan( hypot(xx,yy)/(zz-mean(z)) )/pi*180; % 拟合线和垂线的夹角
    clear direction
end

