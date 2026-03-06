clc;clear;close all;
% 选择样本，读取并绘制原始的SWOT参数与对应的NCOM数据
region   = 'amseas'; % 涡旋所在海域  alaska, amseas
model_path= ['D:\ModelData\NCOM_model\',region,'\'];     % 储存模式数据的路径
eddypath = ['E:\SWOT_NCOM\NCOMGrid\',region,'\'];
outpath  = ['E:\SWOT_NCOM\NCOMEddy\',region,'\'];
% ---------------------------读取NCOM数据列表----------------------------
%{
eddydata = dir([eddypath,'*.mat']); %
for n = 1:length(eddydata)
    load([eddypath,eddydata(n).name]);
    model_file(n,:) = ModelData.NCOM_Name;
end
model_file = unique(model_file,'row');
save(['E:\SWOT_NCOM\NCOMGrid\','ModelFile_',region,'.mat'],'model_file');
clear eddydata n ModelData
%}
load(['E:\SWOT_NCOM\NCOM_match_SWOT\NCOMGrid\','ModelFile_',region,'.mat']);
%% -------------------NCOM三维涡旋检测与垂向追踪--------------------
res = 3.3;
hpkm = 50;  % 对模式数据高通滤波的尺寸
OWValue = -0.2; 
ShapeError = 0.35;% OW参数法检测涡旋的输入
for n=1301:length(model_file)%=1435
% n=1;
    disp(n);
    model_nm = [model_path, model_file(n,:)];
    Model = ReadNCOMmodel([model_path, model_file(n,:)],res,hpkm);% 读model文件
    fields = fieldnames(Model);
    [lat_m, lon_m] = meshgrid(Model.lat, Model.lon);
    for layer=[1,3,6,9,11,13,15,16,18,20:25] % 0,4,10,20,30,40,50,60,80,100,125,150,200,250,300m共15层
%% NCOM Vorticity & Okubo-Weiss 参数计算与轮廓识别
        u_hp = Model.water_u_hp(:,:,layer);
        v_hp = Model.water_v_hp(:,:,layer);
        % ----相对涡度----
        [dudx,dudy] = gradient(permute(u_hp,[2,1,3]),(1/30)/180*pi*6371*10^3); % u
        [dvdx,dvdy] = gradient(permute(v_hp,[2,1,3]),(1/30)/180*pi*6371*10^3); % v
        vorticity = permute(dvdx-dudy,[2,1,3]); % 逆时针为+，顺时针为-
%         zeta = vorticity./gsw_f(lat_m);% 利用f归一化   
        % -----应变率-----
        Sn = permute(dudx-dvdy,[2,1,3]); % 应变的法向分量
        Ss = permute(dvdx+dudy,[2,1,3]); % 应变的剪切分量
        W = Sn.^2 + Ss.^2 - vorticity.^2;% 计算OW参数
        W_norm = zscore_normalize(W); %(基于全场的标准差)归一化涡度的范围选取！！！
        % MSR=permute(hypot(dudx-dvdy,dvdx+dudy),[2,1,3])./gsw_f(lat_m);
        clear u_hp v_hp Sn Ss W dudx dudy dvdx dvdy
        % -------------- 平滑滤波和提取轮廓（阈值0W<-0.2）-----------------
        W_noiseless = filt2(W_norm, 1/30, 1/12, 'lp');
        % 提取涡旋极性（涡度）、位置、轮廓、拟合半径、形状误差
        eddy_contour = get_ow_eddy(lon_m, lat_m, W_noiseless, vorticity, OWValue, ShapeError); 
        depstr = ['Depth_',num2str(Model.depth(layer)),'m'];
        clear ii W W_norm 
        % -------------------- 提取涡旋对应的Model参数 --------------------
        for i=1:length(eddy_contour)
            Cutlim = [min(eddy_contour(i).edgelon), max(eddy_contour(i).edgelon); ...
                min(eddy_contour(i).edgelat), max(eddy_contour(i).edgelat)]; % SWOT Eddy空间范围
            % 模式数据与Eddy匹配，裁剪
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
            eddy_contour(i).lon = Model.lon(lon_i(1):lon_i(2)); % 经度
            eddy_contour(i).lat = Model.lat(lat_i(1):lat_i(2)); % 纬度
            [lat_xy,~] = meshgrid(eddy_contour(i).lat,eddy_contour(i).lon);
            for ff=5:length(fields) % 10为垂向速度
                if layer <= size(Model.(fields{ff}),3)
                    eddy_contour(i).(fields{ff}) = Model.(fields{ff})(lon_i(1):lon_i(2),lat_i(1):lat_i(2),layer);
                end
            end
            eddy_contour(i).vorticity = vorticity(lon_i(1):lon_i(2),lat_i(1):lat_i(2));
            eddy_contour(i).ow = W_noiseless(lon_i(1):lon_i(2),lat_i(1):lat_i(2));
            clear Cutlim lon_i lat_i lat_xy ff dudx dudy dvdx dvdy
        end
        EddyData.(depstr) = eddy_contour;
    end
%% 三维涡旋垂向追踪（判断不同层级的涡旋是否为同一个）
    % 流程：1）极性相同 2）检索范围：中心距离<(R1+R2)  3）若有多个：距离更接近的！
    field = fieldnames(EddyData);
    for ii = 1:length(EddyData.(field{1})) % 单个涡旋
        EddyData.(field{1})(ii).track = ii;
    end
    num = ii; % 第一层最后一个track number
    clear ii
    for dd = 2:length(field) % 层级
        for ii = 1:length(EddyData.(field{dd})) % 单个涡旋
            eddy = EddyData.(field{dd})(ii); % 提取单个涡旋参数
            % 上一层的参数
            type = [EddyData.(field{dd-1}).type]';
            center = cell2mat({EddyData.(field{dd-1}).center}');
            radius = [EddyData.(field{dd-1}).radius]';
            % 匹配（提取上层涡旋参数）
            type_index = find(type==eddy.type); % 极性匹配
            dis = distance(eddy.center(2), eddy.center(1), center(:,2), center(:,1))./180.*pi.*6371; %单位：km 
            rr_index = find(dis-radius<eddy.radius); % 距离匹配
            index = intersect(type_index,rr_index); % 满足1）极性相同 2）检索范围：中心距离<(R1+R2)
            if isempty(index) % 若没有搜寻到上层匹配的，则 num+1 ！
                num = num+1;
                EddyData.(field{dd})(ii).track = num; % 编号
            else
                if length(index) > 1
                    [~,nearest]=min(dis(index)); %3）若有多个：距离最近的！
                    index = index(nearest);
                end
                EddyData.(field{dd})(ii).track = EddyData.(field{dd-1})(index).track;
            end
            clear eddy type center radius type_index dis rr_index index nearest
        end
    end
    save([outpath, '\Eddy_', model_file(n,1:end-3), '.mat'], 'EddyData');
end