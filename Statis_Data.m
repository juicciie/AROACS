function [Best_val,Worst_val,Mean_val,Median_val,Std_val,Wilcoxon_RankSum_Test,Wilcoxon_SignedRank_Test,Friedamn_val] = Statis_Data(result,Func_name,num)

len = length(result.name);
[k,maxIter] = size(result.data{1});
data = zeros(len,k); %保存数据
for i = 1:len
    data(i,:) = min(result.data{i}, [], 2)';
end

Best_val = zeros(1,len); %最优值
Worst_val = zeros(1,len); %最差值
Mean_val = zeros(1,len); %均值
Median_val = zeros(1,len); %中位数
Std_val = zeros(1,len); %标准差

for i = 1:len
    Best_val(i) = min(data(i,:));
    Worst_val(i) = max(data(i,:));
    Mean_val(i) = mean(data(i,:));
    Median_val(i) = median(data(i,:));
    Std_val(i) = std(data(i,:));
end
%-----Wilcoxon秩和检验
Wilcoxon_RankSum_Test = [];
for i = 1:len
    if i == num
        continue;
    end
    [p,h] = ranksum(data(num,:),data(i,:));
    Wilcoxon_RankSum_Test = [Wilcoxon_RankSum_Test;[p h]];
end
%-----Wilcoxon符号秩检验
Wilcoxon_SignedRank_Test = [];
for i = 1:len
    if i == num
        continue;
    end
    [p,h] = signrank(data(num,:),data(i,:));
    Wilcoxon_SignedRank_Test = [Wilcoxon_SignedRank_Test;[p h]];
end
%-----Friedman检验
[p,tbl,stats] = friedman(data');
Friedamn_val = stats.meanranks;
end