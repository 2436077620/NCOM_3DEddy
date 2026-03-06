clc;clear;close all
% 时间序列变化
% 分析正常涡/异常涡的数量、半径、能量、深度、海面风等特征随时间的变化
% 将提取出的NCOM涡旋分为四类：CN AN CA AA
region  = 'alaska'; % 涡旋所在海域  alaska, amseas, us_east
outpath = 'E:\Figure\';
load(['E:\NCOMEddy\','Layer_Number_',region,'.mat']); % 各层级的数量数据
load(['E:\NCOMEddy\','Track_para_',region,'.mat']); % 读取参数的统计数据
PDA(PDA>1 | PDA<-1)=nan;
radius(radius==0)=nan;
dep_even = -300:0;
%% 正常涡/异常涡的数量
% Index
Index.C_norm   = find(type=='C' & flag==1);
Index.A_norm   = find(type=='A' & flag==1);
Index.C_abnorm = find(type=='C' & flag==0);
Index.A_abnorm = find(type=='A' & flag==0);
field = fieldnames(Index);
%% 数量统计随时间的变化
for i = 1:4
    Index.percent(i) = length(Index.(field{i}))/length(type)*100; % 表格
end
% 随月份的占比变化
for m=1:12
    Index.percent_m(m,1) = sum(type=='C' & flag==1 & month==m)/sum(month==m)*100;
    Index.percent_m(m,2) = sum(type=='A' & flag==1 & month==m)/sum(month==m)*100;
    Index.percent_m(m,3) = sum(type=='C' & flag==0 & month==m)/sum(month==m)*100;
    Index.percent_m(m,4) = sum(type=='A' & flag==0 & month==m)/sum(month==m)*100;
end
% ----------
figure;
set(gcf,'position',[200,100,1000,600],'color','w');
ggplotAxes2D(gca,'AxesTheme','ownl','LegendStyle','ownl');
set(gca,'FontName','Times New Roman','FontSize',24,'FontWeight','bold');hold on;
p1 = plot(1:12, Index.percent_m(:,1),'-','linewi',1.8,'color','b');
p2 = plot(1:12, Index.percent_m(:,2),'-','linewi',1.8,'color','r');
p3 = plot(1:12, Index.percent_m(:,3),'--','linewi',1.8,'color','b');
p4 = plot(1:12, Index.percent_m(:,4),'--','linewi',1.8,'color','r');
legend([p1,p2,p3,p4],{'NCE','NAE','ACE','AAE'},'FontSize',16,'location','southwest','color','w');
xlabel('Month');
ylabel('Percentage (%)');
set(gca, 'box', 'on','TickDir','in','xlim',[1 12],'xtick',1:12,'ylim',[0 45],'color','w')
% print(gcf,'-djpeg','-r450',[outpath,'Percentage_month_',region,'.jpg']);