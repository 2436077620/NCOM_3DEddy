clc;clear;close all
% 数量时间变化
% 不同海域的数量随月份变化
%% ---------2024年：涡旋的 数量 随month\depth的变化（柱状图）---------
load(['E:\SWOT_NCOM\NCOMEddy\Track_ECD_Offset\','Track_para_alaska.mat']);
month_num1 = tabulate(month(85963:488764)); % alaska
month([1:85963,488764:end]) = nan;
MLD_mm=[];
for m=1:12
    MLD_mm(m,:) = nanmean(MLD(month == m,:)); % MLD随月份
end
mld1 = mean(MLD_mm,'all','omitnan');

load(['E:\SWOT_NCOM\NCOMEddy\Track_ECD_Offset\','Track_para_amseas.mat']);
month_num2 = tabulate(month(123893:513162)); % amseas
month([1:123893,513162:end]) = nan;
MLD_mm=[];
for m=1:12
    MLD_mm(m,:) = nanmean(MLD(month == m,:)); % MLD随月份
end
mld2 = mean(MLD_mm,'all','omitnan');
clear month MLD_mm
month_num1(:,2) = month_num1(:,2)/10000;
month_num2(:,2) = month_num2(:,2)/10000;
%% Month
figure;
set(gcf,'position',[200,100,1200,450],'color','w');
ggplotAxes2D(gca,'AxesTheme','ownl','LegendStyle','ownl');
set(gca,'FontName','Times New Roman','FontSize',22,'FontWeight','bold');hold on;
box on; grid on;
% 创建柱状图
h2 = bar(month_num2(:,1)-0.18,month_num2(:,2),0.36,'EdgeColor','none','FaceColor',[.1 .5 .5]);
h1 = bar(month_num1(:,1)+0.18,month_num1(:,2),0.36,'EdgeColor','none','FaceColor',[.3 .3 .3]);
alpha(0.8) % 透明度
xlabel('Month');
ylabel('Number (\times10^4)');
legend([h2,h1],{'GoM and Caribbean Sea','Northeast Pacific'},'FontSize',16);
set(gca,'TickDir','in','xlim',[0.4 12.6],'xtick',1:12,'color','w')
print(gcf,'-djpeg','-r450',[outpath,'Number_month.jpg']);