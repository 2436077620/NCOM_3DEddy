function var_stru = NCOM_Eddydec(lon,lat,watervar,varargin)
% 判断NCOM模式是否存在闭合轮廓(基于高通)，用于主函数(NCOMEddyDect_main.m)
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
defaultRadiusThreshold=5; % 默认直径，unit: Km
p.addOptional('radiusThreshold',defaultRadiusThreshold,...
    @(x)validateattributes(x,{'numeric'},...
    {'nonzero','positive'},'eddyDetection_output','radiusThreshold',4));
defaultShapeError=80; % 默认形状误差
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
ShapeError=p.Results.ShapeError; % 形状误差
ARC=6371.393; % 赤道半径
EddyLatDis = (radThreshold/(ARC.*pi/180));                                % 纬度差 直接转为弧度
EddyLonDis = (abs(radThreshold./(ARC.*cos(mean(Y,'all').*pi/180).*pi/180))); % 经度差 与所在纬度有关
layer = [1,6,9,11,13,15:21]; % depth = [0:10:100] 对应 layer = [1,6,9,11,13,15:20]
var_stru = struct();
for dep=1:length(layer) % a 3-D matrix. 
    var_2=watervar(:,:,layer(dep))*100; % sea level anomaly unit from meter to centimeter
    cf = figure('visible','off');
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
                clear dis
                % 基于顶点到圆心距离的偏差计算
                dis=distance(y,x,mean(y),mean(x))/180*pi*6371; %km
                radius = mean(dis);
                % 判据
                E1=abs(dis-radius)>ShapeError * radius; %查找是否存在顶点距离误差>70%半径的位置（1则排除）
                if sum(E1)==0 
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
            [extre_level,ai,~] = unique(kex(:,1:4),'rows');
            var_stru(dep).extre=[extre_level,kex(ai,5)]; % 提取出所有轮廓内的极值（筛除被包含的轮廓）
            % 提取轮廓
            k = k(rowsindex);
            var_stru(dep).edge=k(ai);
        end
    end
    close all;
end
end
%{
% Test Fig
figure;
set(gca,'FontName','Times New Roman','FontSize',24,'FontWeight','bold');hold on;
i=23;
depth = [0:10:100,125,150,200];

[X,Y,Z]=meshgrid(ModelData.Grid(i).lon, ModelData.Grid(i).lat, -ModelData.Grid(i).depth);
xslice=[];
yslice=[];
zslice=-depth(1:2:11);

cen1 = ceil(size(X,1)/2);
cen2 = ceil(size(X,2)/2);
plot3(repmat(X(cen1,cen2,1),1,41), repmat(Y(cen1,cen2,1),1,41), [20;squeeze(Z(cen1,cen2,:))], '-.','linewi',1.8,'color',[.5 .5 .5]);
plot3(repmat(X(cen1,cen2,1),1,41), repmat(Y(cen1,cen2,1),1,41), [20;squeeze(Z(cen1,cen2,:))], '.','Markersize',15,'color',[.5 .5 .5]);

% ----------SWOT------------
height = 20;
[lat_re,lon_re] = meshgrid(SWOTData.Grid(i).lat,SWOTData.Grid(i).lon);
surf(lon_re,lat_re,height.*ones(size(lat_re)),SWOTData.Grid(i).ssha_sub); shading interp% 在z=h的位置绘制图层
colormap(m_colmap('diverging'));  
% caxis([-40 -10]);    
freezeColors;
% 水平流场
val = 2; % 绘制流的像素间隔
u = SWOTData.Grid(i).ua(1:val:end,1:val:end);
v = SWOTData.Grid(i).va(1:val:end,1:val:end);
w = zeros(size(u));  
quiver3(lon_re(1:val:end,1:val:end),lat_re(1:val:end,1:val:end),height.*ones(size(u)),...
 u, v, w, 1,'color',[.6 .6 .6],'linewidth',1.2);
edgezz = ones(size(SWOTData.Eddy(i).edge(1,:)));
plot3(SWOTData.Eddy(i).edge(1,:), SWOTData.Eddy(i).edge(2,:), height.*edgezz,'-b');

% ----------NCOM------------
var = permute(ModelData.Grid(i).temp_hp,[2,1,3]);
var_stru = modelEddy.var_stru_1;

slice(X,Y,Z,var,xslice,yslice,zslice);shading interp
colormap(m_colmap('diverging'));
h1=colorbar;
set(get(h1,'label'),'string','Temperature_h_p (\circC)','FontName','Times New Roman','FontSize',30,'FontWeight','bold'); %Vorticity
caxis([-0.6 0.6]);

% cs = contourslice(X,Y,Z,var,xslice,yslice,zslice,5);
% cs = contourslice_linecolor(cs,'--',0.5,[.6 .6 .6],0.8); %基于contourslice，修改linestyle,linewidth,linecolor,linealpha

for ddd= 1:2:11%length(modelEddy.var_stru_1)
    for nnn = 1:size(var_stru(ddd).edge,2)
        xx = var_stru(ddd).edge(nnn).x;
        yy = var_stru(ddd).edge(nnn).y;
        zz = -depth(ddd).*ones(size(xx));
        plot3(xx, yy, zz,'-k');
    end
end
set(gca,'xlim',[X(1),X(end)],'ylim',[Y(1),Y(end)],'zlim',[-101,20.1],'ztick',-100:20:0,'color','w')
view(-22,10);
xlabel('Longitude');
ylabel('Latitude');
zlabel('Depth(m)');

%}
