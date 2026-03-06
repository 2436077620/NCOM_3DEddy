clc;clear;close all
% （子区域）正常气旋/异常气旋/正常反气旋/异常反气旋的3D复合结构——画图
% 绘制（三维结构的）涡度 & 流场
region = 'amseas'; % 涡旋所在海域  alaska, amseas, us_east
load('E:\DATA\NCOMEddy\','Layer_Number_'); % 各层级的数量数据+各数据的时间数据 
clear Number 
% 复合涡数据路径
path0 = 'E:\SWOT_NCOM\NCOMEddy\Composite_new\';
outpath = 'E:\Figure\';
% 复合涡的输出变量
var_name  = {'water_u_hp';'water_v_hp';'potential_density_hp';'vorticity'};
eddyfield = {'NCE';'NAE';'ACE';'AAE'};
out_name  = {'Velocity';'PDA';'Vorticity';'EKE'};
load([path0,'Composite_2024_',region,'.mat']);
ddep = [1,9,13,16,18,20,21:23];
zbottom = -200;
posit = [650, 0.53, 0.77]; % gcf宽度，gca右边界，colorbar左边界
dep_even = -(0:300);
Depth = -[0:2:12,15:5:50,60:10:100,125,150:50:300];
%% ------------ 3D Figures(全部时间：四个类型) ---------------
for sub = 1:length(Composite) % 不同子区域
    CompEddy = Composite(sub); % 子区域的涡旋
    for ee = 1:4 % 不同涡旋类型
        Num = CompEddy.(eddyfield{ee}).number;
        UA  = CompEddy.(eddyfield{ee}).water_u_hp./Num;
        VA  = CompEddy.(eddyfield{ee}).water_v_hp./Num;
        Vorticity = CompEddy.(eddyfield{ee}).zeta./Num;               % 涡度
        vort_x = permute(Vorticity(41,:,1:25),[3,2,1]);
        vort_y = permute(Vorticity(:,41,1:25),[3,1,2]);
        % ------------------------ 画图 -------------------------
        figure ('visible','off');
        set(gcf,'position',[200,50,posit(1),950],'color','w');
        set(gca,'position',[0.21,0.1,posit(2),0.88],'FontName','Times New Roman','FontSize',26,'FontWeight','bold');
        hold on;
        for dep = ddep
            surf(xx,yy,Depth(dep).*ones(size(xx)),Vorticity(:,:,dep)); shading interp
        end
        % 绘制流速剪头
        res = 1:8:81;
        [XX,YY,ZZ] = meshgrid(xx(1,res), yy(res,1), Depth(ddep));
        quiver3(XX, YY, ZZ, UA(res,res,ddep), VA(res,res,ddep),...
            zeros(length(res),length(res),length(ddep)),0.15,'color',[.4 .4 .4],'linewidth',1);
        colormap(flip(othercolor('Spectral10')));
        clabel = '\zeta/f';
        if ee==1 || ee==3
            caxis([0 0.8])
        else
            caxis([-0.6 0])
        end
        c = colorbar('location','east','position',[posit(3) 0.142 0.034 0.7],'FontSize',24);
        c.Label.String = clabel;
        view(-21,6);
        set(gca,'xlim',[-1.6 1.6],'ylim',[-1.6 1.6],'zlim',[zbottom-0.1,0.1],'ztick',zbottom:50:0,'color','w');
        pos = axis;
        xlabel('\DeltaX')
        ylabel('\DeltaY')
        zlabel('Depth (m)');
        box on
        grid on
        nm = [upper(region(1)),region(2:end)];
        print(gcf,'-djpeg','-r450',[outpath,'3D_vort\Sub',num2str(sub),'_',nm,'_',eddyfield{ee},'_3D.jpg']);
        close all
    end
end