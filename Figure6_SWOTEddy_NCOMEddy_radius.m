clc;clear;close all
% SWOT Eddy和NCOM Eddy时空匹配结果分析
% 用于Model validation
region    = 'alaska'; % 涡旋所在海域  alaska, amseas
load(['E:\DATA\Eddy_OW_Dataset\Vertical_para_',region,'.mat']);
radius(radius>40 | radius<1) = nan;
r_swot = RR;
rr = [];
for n = 1:size(RR,1)
    [~,index] = min(abs(radius(n,:) - RR(n)));
    rr(n,1) = radius(n,index);
end
r_ncom = rr;
%%
region    = 'amseas'; % 涡旋所在海域  alaska, amseas
load(['E:\DATA\Eddy_OW_Dataset\Vertical_para_',region,'.mat']);
radius(radius>40 | radius<1) = nan;
r_swot = [r_swot; RR];
rr = [];
for n = 1:size(RR,1)
    [~,index] = min(abs(radius(n,:) - RR(n)));
    rr(n,1) = radius(n,index);
end
r_ncom = [r_ncom; rr];
% ------------------比较匹配涡旋的半径和振幅-------------------
CData = density2C(r_swot, r_ncom, 0:0.5:20, 0:0.5:20, flip(othercolor('Spectral10'))); close all;
Bias  = mean(r_ncom - r_swot); % 偏差
RMSE = sqrt(mean((r_ncom - r_swot).^2));% 均方根误差
MAPE = mean(abs((r_ncom - r_swot)./r_ncom))*100; % 平均绝对百分比误差
R = corr(r_swot, r_ncom);
N = length(r_swot);
%%
figure;
set(gcf,'position',[200,50,800,650],'color','w');
set(gca,'FontName','Times New Roman','FontSize',22,'FontWeight','bold');hold on;
scatter(r_swot, r_ncom,'filled','CData',CData);
colormap(flip(othercolor('Spectral10')));
colorbar;
plot([0 100], [0 100], '-k','linewi',1);
text(1,19,['R = ',num2str(R,'%.2f')],'FontName','Times New Roman','FontSize',18,'FontWeight','bold');
text(1,17.5,['RMSE = ',num2str(RMSE,'%.2f'),' km'],'FontName','Times New Roman','FontSize',18,'FontWeight','bold');
text(1,16,['Bias = ',num2str(Bias,'%.2f'),' km'],'FontName','Times New Roman','FontSize',18,'FontWeight','bold');
text(1,14.5,['MAPE = ',num2str(MAPE,'%.2f'),' %'],'FontName','Times New Roman','FontSize',18,'FontWeight','bold');
text(1,13,['N = ',num2str(N)],'FontName','Times New Roman','FontSize',18,'FontWeight','bold');
xlabel('Radius_{\fontsize{16}SWOT eddy} (km)');
ylabel('Radius_{\fontsize{16}NCOM eddy} (km)');
xlim([0 20]);ylim([0 20]);
yticks(0:5:20);
box on; grid on;
% print(gcf,'-djpeg','-r550',['E:\Figure\Radius2.jpg']);
