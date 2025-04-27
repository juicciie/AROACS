clear,clc,close all
format long
N = 50; 
Max_iter = 1000; 
runs = 30; 
% EFO default configurations
c1 = 0.26;
c2 = 1.26;
Pout = 0.31;
result = {}; 
name_list = {'AROACS','AROA','AEFA','ASO','EFO','GSA','KOA','POA','RPSO','WSAR'}; 
len = length(name_list);
result.name = name_list; 
result.data = cell(1,len); 
friedman_path = '/Users/juicciie/Paper/CEC 2017/original data/100D/friedman/';
wilcoxon_ranksum_path = '/Users/juicciie/Paper/CEC 2017/original data/100D/wilcoxon_ranksum/';
wilcoxon_signedrank_path = '/Users/juicciie/Paper/CEC 2017/original data/100D/wilcoxon_signedrank/';
for i = 1:len
    result.data{i} = zeros(runs,Max_iter);
end
for k = 3:30
    Function_name = sprintf(['F', num2str(k)]);
    [lb, ub, dim, fobj] = CEC2017_100(Function_name);
    for i = 1:runs
        [Best_fitness1, Best_pos1, Iter_curve1] = AROA_CS(N, Max_iter, lb, ub, dim, fobj);
        [Best_fitness2, Best_pos2, Iter_curve2] = AROA(N, Max_iter, lb, ub, dim, fobj);
        [Best_fitness3, Best_pos3, Iter_curve3] = AEFA(N, Max_iter, lb, ub, dim, fobj, 1, 1, 2);
        [Best_pos4, Best_fitness4, Iter_curve4] = ASO(N, Max_iter, lb, ub, dim, fobj);
        [Best_pos5, Best_fitness5, Iter_curve5] = EFO(fobj, dim, lb, ub, N, Max_iter, c1, c2, Pout);
        [Best_fitness6, Best_pos6, Iter_curve6] = GSA(N, Max_iter, lb, ub, dim, fobj, 1, 1, 2);
        [Best_fitness7, Best_pos7, Iter_curve7] = KOA(N, Max_iter, lb, ub, dim, fobj);
        [Best_pos8, Best_fitness8, Iter_curve8] = POA(N, Max_iter, lb, ub, dim, fobj);
        [Best_pos9, Best_fitness9, Iter_curve9] = RPSO(N, Max_iter, lb, ub, dim, fobj);
        [Best_fitness10, Best_pos10, Iter_curve10] = WSAR(N, Max_iter, lb, ub, dim, fobj);
    
        % 存储运算结果
        result.data{1}(i,:) = Iter_curve1;
        result.data{2}(i,:) = Iter_curve2;
        result.data{3}(i,:) = Iter_curve3;
        result.data{4}(i,:) = Iter_curve4;
        result.data{5}(i,:) = Iter_curve5;
        result.data{6}(i,:) = Iter_curve6;
        result.data{7}(i,:) = Iter_curve7;
        result.data{8}(i,:) = Iter_curve8;
        result.data{9}(i,:) = Iter_curve9;
        result.data{10}(i,:) = Iter_curve10;
    end
    num = 1; %选择算法num和其他算法进行对比（秩和检验，符号秩检验）
    [Best_val,Worst_val,Mean_val,Median_val,Std_val,Wilcoxon_RankSum_Test,...
        Wilcoxon_SignedRank_Test,Friedamn_val] = Statis_Data(result,Function_name,num);
    save([friedman_path, Function_name], 'Friedamn_val');
    save([wilcoxon_ranksum_path, Function_name], 'Wilcoxon_RankSum_Test');
    save([wilcoxon_signedrank_path, Function_name], 'Wilcoxon_SignedRank_Test');
end
