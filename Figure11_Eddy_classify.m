clc;clear;close all
% 基于OW检测的 NCOMData 计算涡核深度
% UA、VA即为涡致水平流速
region    = 'amseas'; % 涡旋所在海域  alaska, amseas, us_east
ncom_eddy = 'E:\DATA\NCOMEddy\'; 
model_path= 'E:\DATA\NCOM_model\';     % 储存模式数据的路径
outpath = 'E:\Figure\';
model_file = dir([ncom_eddy,'Eddy*.mat']);
res = 3.3;
hpkm = 50;  % 对模式数据高通滤波的尺寸
%% 读取全局模式参数 & 检测出的NCOM涡旋
n = 1;
    Model = ReadNCOMmodel([model_path, model_file(n).name(6:end-4),'.nc'],res,hpkm);% 读model文件
    load([ncom_eddy,model_file(n).name]); % 读取对应的NCOMEddy
    field = fieldnames(EddyData);
    [Count,Depth] = NCOM_Eddytrack(EddyData,3);
%% 从Model中提取信号（与涡旋信号比较），lonlim & latlim的范围    
num = 110; % 样本——n=1，num=158:反气旋; 26:气旋; 38（183）:反气旋双峰; 189（110）:异常气旋 48:异常反气旋
    if num == 158
        name = 'normA';
    elseif num == 26
        name = 'normC';
    elseif num == 38
        name = 'doubleA';
    elseif num == 110
        name = 'abnormC';
    elseif num == 48
        name = 'abnormA';
    end
    % ---------------------- 提取涡旋对应的Model参数 ----------------------
    lonlim = []; latlim = []; % 提取裁剪范围
    for ii = 1:size(Count(num).layer,1)
        index = Count(num).layer(ii,2);
        eddy_para = EddyData.(field{Count(num).layer(ii,1)})(index);
        lonlim = [min([eddy_para.edgelon, lonlim]),max([eddy_para.edgelon, lonlim])]; % 范围
        latlim = [min([eddy_para.edgelat, latlim]),max([eddy_para.edgelat, latlim])]; % 范围
    end
    Cutlim = [lonlim(1), lonlim(2); latlim(1), latlim(2)]; % SWOT Eddy空间范围
    Eddymodel = cutNCOMgrid(Model,Cutlim); % cutNCOMgrid用于裁剪和提取参数
    clear Cutlim index lonlim latlim ii eddy_para
    % ----------计算PDA（基于DA变化的统计关系，需要Model参数）-------------
    PDA = NCOM_EddyPDA(EddyData, Eddymodel, Count(num));
    if num == 110
        PDA(PDA>0) = PDA(PDA>0)./8;
    end
    %% --------------------平滑PDA & 计算涡核深度ECD------------------------
    dep_even = Depth(1):-1:Depth(end); 
    x1 = interp1(Depth, PDA, dep_even, 'pchip'); % 插值到1m间隔
    x1 = smooth(x1, 31);       % 21点平滑滤波器（相当于考虑20m的深度层）
    xx = zscore(x1); % PDA归一化
    for fl = 1:1
        if strcmp(name,'normA') || strcmp(name,'doubleA') % 理论上：反气旋存在PDA极小值
            [a,dep] = findpeaks(-xx); % 查找极小值
            Npda = -x1(dep);
            xrange = [-1 0 0 -1];
        elseif strcmp(name,'normC')
            [a,dep] = findpeaks(xx); % 查找极大值
            Npda = x1(dep);
            xrange = [1 0 0 1];
        elseif strcmp(name,'abnormC')
            [a,dep] = findpeaks(-xx); % 查找极小值
            Npda = x1(dep);
            xrange = [-1 0 0 -1];
        elseif strcmp(name,'abnormA')
            [a,dep] = findpeaks(xx); % 查找极大值
            Npda = -x1(dep);
            xrange = [1 0 0 1];
        end
    end
    % ----------------------计算混合层深度MLD-------------------------
    T = [];
    for dd = 1:25
        T(dd,1) =  mean(Eddymodel.temp(:,:,dd),'all','omitnan'); % 网格平均温度
    end
    dep_5 = 0:-1:-300; 
    tt = interp1(-Eddymodel.depth(1:25), T, dep_5, 'pichp'); % 插值到1m间隔
    tt_index = find(tt<T(6)-0.2);
    MLD_dep = dep_5(tt_index(1));
    clear dd T dep_5 tt tt_index
%% 画图
figure;
if Count(num).type == 'A'
    color = [1 0 0];
    face = [.9 .5 .5];
    location = 'southwest';
else
    color = [0 0 1];
    face = [.5 .5 .9];
    location = 'southwest';
end
set(gcf,'position',[200,50,600,900],'color','w');
set(gca,'fontsize',24,'Fontname','Time new roman','fontweight','bold');hold on;
xline(0,'-k','linewi',1.2); % PDA的0线
% 基于OW参数法识别-垂向追踪的三维涡旋 
top=Count(num).layer(1,1);
bottom=Count(num).layer(end,1);
yline(Depth(top),'-','linewi',1.5,'color',color);
yline(Depth(bottom),'-','linewi',1.5,'color',color); % 涡旋深度（最上层-底层）
L1 =fill(xrange,[Depth(top),Depth(top),Depth(bottom),Depth(bottom)],face,'FaceAlpha',0.2,'EdgeAlpha',0);
L2 = yline(MLD_dep,'-.','linewi',2,'color',[.2 .7 .2]); % MLD
% 三维涡旋平均垂向PDA
p1 = plot(PDA,Depth,'--','linewi',1.8,'color',[.4 .4 .4]);
p2 = plot(x1,dep_even,'-','linewi',2,'color',color);
p3 = plot(Npda,dep_even(dep),'k.','Markersize',30,'linewi',2.5);
xlabel('PDA(kg/m^3)');
ylabel('Depth (m)');
xlim([min(PDA)-0.0065,max(PDA)+0.0075]);
ylim([-300 0]);  
legend([p1,p2,L2,L1],{'Original signal','Filtered signal','MLD','OW eddy range'},'location',location);
box on
grid on
% print(gcf,'-djpeg','-r450',[outpath,'Example_EddyCore_',name,'2.jpg']);