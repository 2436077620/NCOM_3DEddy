function var_stru = NCOM_Surfel_dec(lon,lat,watervar,varargin)
% 判断NCOM模式的surf_el_hp是否存在闭合轮廓(基于高通)，用于主函数(NCOM_Eddy_Dataself_self.m)
    warning('off')			 %根据这个唯一标识符隐藏这条warning
    % % Inputs validation process
    if nargin<3
        error('Not Enough Inputs. 3 inputs are required.');
    end
    if nargin>5
        error('Too Many Inputs. No more than 4 inputs please.');
    end

    p = inputParser;
    p.addRequired('lon',@(x)validateattributes(x,{'numeric'},...
        {'nonempty'},'eddyDetection_output','lon',1));
    p.addRequired('lat',@(x)validateattributes(x,{'numeric'},...
        {'nonempty'},'eddyDetection_output','lat',2));
    p.addRequired('watervar',@(x)validateattributes(x,{'numeric'},...
        {'3d'},'eddyDetection_output','watervar',3));
    % -----------------------------varargin------------------------------------
    defaultRadiusThreshold=3; % 默认半径，unit: Km
    p.addOptional('radiusThreshold',defaultRadiusThreshold,...
        @(x)validateattributes(x,{'numeric'},...
        {'nonzero','positive'},'eddyDetection_output','radiusThreshold',4));
    defaultShapeError=55; % 默认形状误差
    p.addOptional('ShapeError',defaultShapeError,...
        @(x)validateattributes(x,{'numeric'},...
        {'nonzero','positive'},'eddyDetection_output','ShapeError',5));

    p.parse(lon,lat,watervar,varargin{:});

    % Test whether "lon" and "lat" are vectors or matrixes.
    if size(p.Results.lon,2)==1 && size(p.Results.lat,2)==1
        [X,Y]=meshgrid(lon,lat);
    elseif size(p.Results.lon,1)==1 && size(p.Results.lat,1)==1
        [X,Y]=meshgrid(lon,lat);
    else
        X=lon;Y=lat;
    end

    % test whether sla marix should be rotated
    if isequal(size(p.Results.watervar(:,:,1)),size(X))==0
        watervar=p.Results.watervar;
        watervar=permute(watervar,[2,1,3]);
    else
        watervar=p.Results.watervar;
    end
    if isequal(size(watervar(:,:,1)),size(X))==0
        error('The zonal or meridional grid points of "sla" do not match those of "lon" or "lat".');
    end

    % specify the area for M_Map toolbox
    m_proj('miller','lon',[min(p.Results.lon(:)) max(p.Results.lon(:))],...
        'lat',[min(p.Results.lat(:)) max(p.Results.lat(:))]);

    radThreshold=p.Results.radiusThreshold; %输入半径
    Error=p.Results.ShapeError; % 形状误差
    ARC=6371.393; % 赤道半径
    EddyLatDis = (radThreshold/(ARC.*pi/180));                                % 纬度差 直接转为弧度
    EddyLonDis = (abs(radThreshold./(ARC.*cos(mean(Y,'all').*pi/180).*pi/180))); % 经度差 与所在纬度有关
    var_2=watervar(:,:)*100; % sea level anomaly unit from meter to centimeter
    var_stru =struct();
    
    figure('visible','off');
    [c,~]=contour(X,Y,var_2,100);  
    s = getcontourlines(c); % getcontourlines 解析contour得到的"闭合"轮廓
    if ~isempty(fieldnames(s))
        % 保留尺寸大于EddyLonDis的闭合轮廓(约为2.5个网格大小)
        ss = 0;
        clear k
        for i=1:length(s) 
            x = sort(s(i).x);   %经度
            y = sort(s(i).y);   %纬度
            if abs(x(1)-x(length(x)))>EddyLonDis && abs(y(1)-y(length(y)))>EddyLatDis
                center = [mean(s(i).x),mean(s(i).y)]; % 拟合圆心
                [arclen,az]=distance(s(i).y,s(i).x,center(2),center(1));
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
                    ss=ss+1;
                    k(ss) = s(i);
                end
            end
        end
        % 提取符合要求的闭合轮廓内极值
        if exist('k','var')
            for j=1:length(k)
                in=inpolygon(X,Y,k(j).x,k(j).y);
                XX=X(in); 
                YY=Y(in);
                if mean(var_2(in),'omitnan')>k(j).v % 1：正异常
                    [k(j).extremum(3,1),lloc] = max(var_2(in),[],'omitnan');
                    k(j).extremum(4,1) = 1;
                elseif mean(var_2(in),'omitnan')<=k(j).v % -1：负异常
                    [k(j).extremum(3,1),lloc] = min(var_2(in),[],'omitnan');
                    k(j).extremum(4,1) = -1;
                end
                k(j).extremum(1:2,1) = [XX(lloc),YY(lloc)];
                k(j).extremum(5,1) = nanmean(distance(nanmean(k(j).x),nanmean(k(j).y),k(j).x,k(j).y))/180*pi*ARC; % 半径(°)
            end
            % 极值相同的点，仅保留最外部轮廓半径
            [kex,rowsindex] = sortrows([k.extremum]','descend'); % 降序排列，第一个为半径最大值                
            [~,ai,~] = unique(kex(:,1:4),'rows');
            % 提取轮廓
            k = k(rowsindex);
            for ii=1:length(ai)
                var_stru(ii).edgelon = k(ai(ii)).x;
                var_stru(ii).edgelat = k(ai(ii)).y;
                var_stru(ii).edgevalue = k(ai(ii)).v;
                var_stru(ii).amp = k(ai(ii)).extremum(1:3,1);
                var_stru(ii).type = k(ai(ii)).extremum(4,1);
                var_stru(ii).radius = k(ai(ii)).extremum(5,1);
            end
        end
    end
    close all;
end