function [Count,Depth] = NCOM_Eddytrack(EddyData,thick)
% 从EddyData（基于OW参数检测出的NCOM涡旋数据集）中提取追踪涡旋的Count
% 输出Count包含（基础参数）：track layer 
    field = fieldnames(EddyData); 
    Depth = [];
    for dd = 1:length(field) % 层级
        Depth(dd,1) = -str2num(field{dd}(7:end-1)); % 深度
        track0 = [EddyData.(field{dd}).track]';
        len0 = length(track0);
        if dd==1
            track = [track0,repmat(dd,len0,1),(1:len0)'];
        else
            track = [track; [track0,repmat(dd,len0,1),(1:len0)'] ];
        end
    end
    track1 = unique(track(:,1)); 
    count = 0;
    Count = struct();
    for i = 1:length(track1)
        ii = find(track(:,1)==track1(i));
        if length(ii)>thick    % 删掉小于xx层的涡旋数据
            if length(track(ii,2)) == length(unique(track(ii,2))) % 确保一层不会出现两个同编号涡旋
                count = count+1;
                Count(count).track = track1(i);
                Count(count).layer = track(ii,2:3);
            end
        end
    end
    if ~isempty(fieldnames(Count))
        % -------Count记录了track对应的layer & index---------
        for nn=1:length(Count) % Count中的第nn个涡
            clear LL index type kk center radius
            LL = Count(nn).layer(:,1); % 层
            index = Count(nn).layer(:,2); % 索引 
            % NCOM Eddy参数
            for kk = 1:length(LL)
                center(kk,:) = EddyData.(field{LL(kk)})(index(kk)).center;
                radius(kk,:) = EddyData.(field{LL(kk)})(index(kk)).radius;
            end
            Count(nn).type = EddyData.(field{LL(1)})(index(1)).type; % 极性
            Count(nn).center = center; % 中心位置
            Count(nn).radius = radius; % 半径
        end
    end
end

