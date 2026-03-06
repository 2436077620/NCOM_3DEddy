clc;clear;close all
% 数量垂向变化
% 不同层的涡旋数量随月份变化
%% Depth 总数量随深度的变化
load(['E:\DATA\NCOMEddy\','Layer_Number_alaska.mat']); % 各层级的数量数据
depth_num1 = sum(Number);% alaska
load(['E:\DATA\NCOMEddy\','Layer_Number_amseas.mat']); % 各层级的数量数据
depth_num2 = sum(Number);% amseas
% --------------amseas----------------
subplot(1,2,1);
set(gca,'FontName','Times New Roman','FontSize',24,'FontWeight','bold');hold on;
p2 = plot(depth_num2/(10^5),Depth,'-o','linewi',2,'color',[.1 .5 .5]);
m2 = yline(mld2,'--','linewi',2,'color',[.1 .6 .5]);
text(3.8,mld2+7,['MLD = ',num2str(-mld2,'%.2f'),'m'],'FontName','Times New Roman','FontSize',18,'FontWeight','bold');
legend([p2],{'GoM and Caribbean Sea'},'FontSize',16,'location','northOutside','Box','off');
xlabel('Number (\times10^5)');
ylabel('Depth (m)');
set(gca,'TickDir','out','xlim',[2.8 15],'ylim',[-300 0],'color','w');
box on; grid on
% --------------alaska----------------
figure;
set(gcf,'position',[200,100,1000,900],'color','w');
subplot(1,2,2);
set(gca,'FontName','Times New Roman','FontSize',24,'FontWeight','bold');hold on;
p1 = plot(depth_num1/(10^5),Depth,'-o','linewi',2,'color','k');
m1 = yline(mld1,'--','linewi',2,'color',[.3 .3 .3]);
text(9.3,mld1+7,['MLD = ',num2str(-mld1,'%.2f'),'m'],'FontName','Times New Roman','FontSize',18,'FontWeight','bold');
legend([p1],{'Northeast Pacific'},'FontSize',16,'location','northOutside','Box','off');
xlabel('Number (\times10^5)');
ylabel('');
set(gca,'TickDir','out','xlim',[7 12],'ylim',[-300 0],'color','w');
box on; grid on
% print(gcf,'-djpeg','-r450',[outpath,'Number_depth.jpg']);