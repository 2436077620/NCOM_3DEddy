clc;clear;close all
% 涡旋的平均PDA/MLD的垂向结构
region  = 'amseas'; % 涡旋所在海域  alaska, amseas, us_east
disp(region);
load(['E:\DATA\Track_para_',region,'.mat']); % 读取参数的统计数据
dep_ratio = repmat(0:350, length(MLD),1) ./ repmat(MLD,1,351); % 比值：Depth / MLD
clear ii
load(['E:\DATA\PDA_MLD_',region,'_2024.mat']);
%% 分类    
NCE = type=='C' & flag==1 & month>4 & month<11;
NAE = type=='A' & flag==1 & month>4 & month<11;
ACE = type=='C' & flag==0 & month>4 & month<11;
AAE = type=='A' & flag==0 & month>4 & month<11;
pda_nce = smooth(mean(pda_int(NCE,:),'omitnan'),5);
pda_nae = smooth(mean(pda_int(NAE,:),'omitnan'),5);
pda_ace = smooth(mean(pda_int(ACE,:),'omitnan'),5);
pda_aae = smooth(mean(pda_int(AAE,:),'omitnan'),5);
% 不等间距纵坐标
dep_unequal = -(1:length(dep_int)); % 设置等间距数值
ratio = -[1, 1.5, 2,5, 10:10:-dep_int(end)]; % 刻度对应深度
[~,x] = ismember(ratio, dep_int);   % 深度对应数值
%% 画图： pda（深度/MLD）随深度的变化
figure;
set(gcf,'position',[200,100,800,800],'color','w');
set(gca,'FontName','Times New Roman','FontSize',24,'FontWeight','bold');hold on;
p1 = plot(pda_nce,dep_unequal,'-','linewi',2,'color','b');
p2 = plot(pda_nae,dep_unequal,'-','linewi',2,'color','r');
p3 = plot(pda_ace,dep_unequal,'--','linewi',2,'color','b');
p4 = plot(pda_aae,dep_unequal,'--','linewi',2,'color','r');
xline(0,'-','linewi',1,'color','k');
yline(-x(1),'--','linewi',1.5,'color',[.3 .3 .3]);
yline(-x(2),'-.','linewi',1.5,'color',[.3 .3 .3]);
yline(-x(3),'--','linewi',1.5,'color',[.3 .3 .3]);
legend([p1,p2,p3,p4],{'NCE','NAE','ACE','AAE'},'FontSize',20,'Box','off','location','southeast');
xlabel('PDA (kg/m^3)');
ylabel('z (depth/MLD)');
set(gca,'TickDir','out','ylim',[-length(dep_int) -1],'color','w');
if strcmp(region,'alaska')
    xlim([-0.035 0.035]);
    xticks(-0.03:0.01:0.03);
    title('Northeast Pacific','FontName','Times New Roman','FontSize',24,'FontWeight','bold');
elseif strcmp(region,'amseas')
    xlim([-0.05 0.05]);
    xticks(-0.04:0.02:0.04);
    title('GoM and Caribbean Sea','FontName','Times New Roman','FontSize',24,'FontWeight','bold')
end
yticks(flip(-x));
yticklabels(flip(-ratio));
box on; 
% print(gcf,'-djpeg','-r450',[outpath,'PDA_ratio_',region,'_2024.jpg']);