clc;clear;close all
% 各个子区域的3D复合结构——画图（经向纬向的PDA切片图）
region = 'amseas'; % 涡旋所在海域  alaska, amseas
nm = [upper(region(1)),region(2:end)];
path0 = 'E:\DATA\Composite\';
outpath = 'E:\Figure\';
load([path0,'Composite_2024_',region,'.mat']); % 读取复合涡数据
eddyfield = fields(Composite);
% 读取深度数据Depth
grid_path = ['E:\DATA\',region,'\'];
grid_file = dir([grid_path,'*.mat']);
load([grid_path,grid_file(1).name]); 
Depth = Eddymodel(1).depth;
clear Eddymodel SSHA_Eddy grid_path grid_file
posit = [680, 0.6, 0.81]; % 绘制X-Depth垂面： gcf宽度，gca右边界，colorbar左边界
%% 读取各个子区域的数据
for sub = 1:length(Composite) % 不同子区域
CompEddy = Composite(sub); % 子区域的涡旋
out_name1 = {'PDA';'Vorticity';'EKE'};
[x1,dep1] = meshgrid(xx(1,:),-Depth(1:25));  % 300m以上
%% for ee=1:4 % 不同涡旋
    for ee = 1:4
        Num = CompEddy.(eddyfield{ee}).number;
        % ------------------------ X-Depth垂面图 -------------------------
        PDA = CompEddy.(eddyfield{ee}).potential_density_hp./Num*100.*cha; % 位势密度异常
        var_x = permute(PDA(41,:,1:25),[3,2,1]);
        var_y = permute(PDA(41,:,1:25),[3,1,2]);
        % ----------------经向---------------------
        figure('visible','off');
        set(gcf,'position',[200,50,posit(1),550],'color','w');
        set(gca,'position',[0.19,0.15,posit(2),0.8],'FontName','Times New Roman','FontSize',22,'FontWeight','bold');hold on;
        pcolor(x1,dep1,var_x);shading interp
        contour(x1,dep1,var_x,6,'LineColor',[.5 .5 .5]);
        colormap(flip(othercolor('RdBu11')));
        clabel = 'PDA (\times10^-^2 kg/m^3)';
        caxis([-5 5])
        setPivot(0);
        c = colorbar('location','east','position',[posit(3) 0.15 0.025 0.8],'FontSize',22);
        c.Label.String = clabel;
        xlim([-2 2]); ylim([-300 0]);
        xlabel('X(r)')
        ylabel('Depth (m)');
        plot([-2 2 2 -2 -2],[-300 -300 0 0 -300],'-k','linewi',1);
%         print(gcf,'-djpeg','-r450',[outpath,'sub',num2str(sub),'\',nm,'_',eddyfield{ee},'_X.jpg']);
        % ----------------纬向---------------------
        figure('visible','off');
        set(gcf,'position',[200,50,posit(1),550],'color','w');
        set(gca,'position',[0.19,0.15,posit(2),0.8],'FontName','Times New Roman','FontSize',22,'FontWeight','bold');hold on;
        pcolor(x1,dep1,var_y);shading interp
        contour(x1,dep1,var_y,6,'LineColor',[.5 .5 .5]);
        colormap(flip(othercolor('RdBu11')));
        clabel = 'PDA (\times10^-^2 kg/m^3)';
        caxis([-5 5])
        setPivot(0);
        c = colorbar('location','east','position',[posit(3) 0.15 0.025 0.8],'FontSize',22);
        c.Label.String = clabel;
        xlim([-2 2]); ylim([-300 0]);
        xlabel('Y(r)')
        ylabel('Depth (m)');
        plot([-2 2 2 -2 -2],[-300 -300 0 0 -300],'-k','linewi',1);
%         print(gcf,'-djpeg','-r450',[outpath,'sub',num2str(sub),'\',nm,'_',eddyfield{ee},'_Y.jpg']);
        close all
    end
end