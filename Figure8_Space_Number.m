clc;clear;close all
% 涡旋数量的空间分布与高发区
region  = 'amseas'; % 涡旋所在海域  alaska, amseas, us_east
outpath = 'E:\Figure\';
load(['E:\DATA\','Layer_Number_',region,'.mat']); % 各层级的数量数据
load(['E:\DATA\','Track_para_',region,'.mat']); % 读取参数的统计数据
dep_even = -300:0;
%% 空间位置分布（地理统计）
if strcmp(region,'alaska')
    range = [190, 240, 36.45, 62.05]; % 0~360
    lonrange = [range(1)+1 range(2)-5];
    latrange = [range(3)+3 range(4)-1.3];
    range1 = [360-160.5,360-154,53,57];
    range2 = [360-142,360-136,53.5,57.5];
    range3 = [360-159.5,360-153,40.5,45];
    range4 = [0,0,0,0];
    range5 = [0,0,0,0];
elseif strcmp(region,'amseas')
    range = [260, 305.1, 8, 33];
    lonrange = [range(1) range(2)-5];
    latrange = [range(3)+1.2 range(4)-3];
    range1 = [360-77.5,360-72.5,20,24];
    range2 = [360-64,360-61,11.9,16.1];
    range3 = [360-86,360-79,18,21.5];
    range4 = [360-89,360-84.5,26,29.2];% 该海域选择了4个高发区
end
% 涡旋分类
Index.All = flag<999;
Index.NCE = type == 'C' & flag==1;
Index.NAE = type == 'A' & flag==1;
Index.ACE = type == 'C' & flag==0;
Index.AAE = type == 'A' & flag==0;
field = fieldnames(Index);
[lat,lon] = meshgrid(range(3):1:range(4), range(1):1:range(2));
Space_num = [];
Space_ecd = [];
for ff=1:1
    for ii=1:size(lat,1)-1          
        for jj=1:size(lat,2)-1
            index = (center(:,1)>=lon(ii,1))&(center(:,1)<lon(ii+1,1))&(center(:,2)>=lat(1,jj))&(center(:,2)<lat(1,jj+1)&...
                Index.(field{ff}));
            Space_num(ii,jj)=sum(index);
            Space_ecd(ii,jj)=mean(ECD(index),'omitnan');
            Space_mld(ii,jj)=mean(MLD(index),'omitnan');
            Space_ratio(ii,jj)=mean(ECD(index)./MLD(index),'omitnan');
        end
    end
%% Figure——Num
    figure;
    set(gcf,'position',[200,50,1050,800],'color','w');
    set(gca,'fontsize',24,'Fontname','Time new roman','fontweight','bold');hold on;  
    m_proj('Miller','longitudes',lonrange,'latitudes',latrange); 
    % m_plot(eddycen(:,1), eddycen(:,2),'.','color',[250,67,67]./255,'markersize',5);
    m_pcolor(lon(1:end-1,1:end-1),lat(1:end-1,1:end-1), Space_num./1000);shading interp;
    caxis([0 2]);
    colormap(flip(othercolor('Spectral10')));
    c = colorbar('fontsize',24,'Fontname','Time new roman','fontweight','bold');
    set(c, 'Ticks', 0:3);
    c.Label.String = 'Number (\times10^3)';
    %---------------------------------海岸线--------------------------------------
    m_gshhs_i('patch',[0.5 0.5 0.5],'linewidth',1,'edgecolor','none');
    m_grid('xtick',5,'box','fancy','tickdir','out','Fontname','Time new roman','fontweight','bold','LineWidth',0.5,'axis','equal');
    % 涡旋高发区
    m_plot([range1(1),range1(1),range1(2),range1(2),range1(1)],...
        [range1(3),range1(4),range1(4),range1(3),range1(3)],'linewi',1.5,'color','k');
    m_plot([range2(1),range2(1),range2(2),range2(2),range2(1)],...
        [range2(3),range2(4),range2(4),range2(3),range2(3)],'linewi',1.5,'color','k');
    m_plot([range3(1),range3(1),range3(2),range3(2),range3(1)],...
        [range3(3),range3(4),range3(4),range3(3),range3(3)],'linewi',1.5,'color','k');
    m_plot([range4(1),range4(1),range4(2),range4(2),range4(1)],...
        [range4(3),range4(4),range4(4),range4(3),range4(3)],'linewi',1.5,'color','k');
%     print(gcf,'-djpeg','-r450',[outpath,'Space\Number_spacial_',region,'_sub_1.jpg']);
end