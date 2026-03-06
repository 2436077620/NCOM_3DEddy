function norm_mat = zscore_normalize(mat)
% Z-score 标准化 (Standardization)
% 将数据转换为均值为0，标准差为1的分布
    mu = mean(mat(:),'all','omitnan');
    sigma = std(mat(:),1,'all','omitnan');
    norm_mat = (mat - mu) / sigma;
end

