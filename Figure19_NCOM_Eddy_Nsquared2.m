clc;clear;close all;
% 区域的层化强度：浮力频率
region   = 'amseas'; % 涡旋所在海域  alaska, amseas
model_path= ['E:\DATA\NCOMEddy\',region,'\'];     % 储存裁剪后模式数据的路径
load(['E:\DATA\Track_para_',region,'.mat']);
file = dir([model_path,'*.mat']);
Nsquared = [];
for n = 1:length(file)
    load([model_path,file(n).name]);
    if n==1
        Depth = -Eddymodel(1).depth(1:end-1);
    end
    for i = 1:length(Eddymodel)
        data = Eddymodel(i);
        SA = squeeze(mean(data.salinity,[1,2],'omitnan')); % Absolute Salinity [g/kg] = [psu]
        CT = squeeze(mean(data.temp,[1,2],'omitnan'));     % Conservative Temperature [deg C]
        [Y,~,Z] = meshgrid(data.lat,data.lon,data.depth);
        p = squeeze(mean(ZP(Z,Y,0),[1,2]));       % sea pressure [dbar] 水深转为水压
        lat = squeeze(mean(Y,[1,2]));
        [N2, ~] = gsw_Nsquared(SA,CT,p,lat);
        Nsquared = [Nsquared;N2'];
        clear data SA CT Y Z p lat N2
    end
    clear Eddymodel
end
%% N2的时间-深度变化
month1 = month;
if strcmp(region,'alaska')
    month1([1:85963,488764:end]) = nan;
elseif strcmp(region,'amseas')
    month1([1:123893,513162:end]) = nan;
end
ECD_mm=[]; N2_mm=[]; MLD_mm=[];
for m=1:12
    ECD_n(m,:) = nanmean(ECD(month1 == m & flag == 1)); % 正常涡随月份
    ECD_a(m,:) = nanmean(ECD(month1 == m & flag == 0));
    N2_mm(m,:) = mean(Nsquared(month1 == m,1:25),'omitnan'); % Nsquared随月份
end
%% --------------------------------------------------
figure;
set(gcf,'position',[200,100,1000,600],'color','w');
ggplotAxes2D(gca,'AxesTheme','ownl','LegendStyle','ownl');
set(gca,'FontName','Times New Roman','FontSize',24,'FontWeight','bold');hold on;
[dep_xy,mon_xy] = meshgrid(Depth(1:25),1:12);
pcolor(mon_xy, dep_xy, N2_mm.*10^4);shading interp
colormap(flip(othercolor('Spectral10')));
h1 = colorbar;
set(get(h1,'label'),'string','N^2 (10^{-4} s^{-2})','FontName','Times New Roman','FontSize',28,'FontWeight','bold');
caxis([0 4]);
xlabel('Month');
ylabel('Depth (m)');
set(gca, 'box', 'off','TickDir','out','xlim',[1-0.04 12],'xtick',1:12,'ylim',[-300-2 0],'color','w')
print(gcf,'-djpeg','-r450',['E:\Figure\N2_month_',region,'.jpg']);