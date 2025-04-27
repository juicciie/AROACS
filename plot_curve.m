clear
clc
close all
addpath input_data

Max_iter = 1000; 
N = 50;
run_times = 30;

% EFO default configurations
c1 = 0.26;
c2 = 1.26;
Pout = 0.31;

Function_name = 'F1';         
Optimal_results={}; 
[lb, ub, dim, fobj] = CEC2017_50(Function_name);
addpath(genpath('comparison for attraction'))

for run_time = 1:run_times
% -----------------------------------你的算法放在首位--------------------------------
tic
[Best_f,Best_x,cg_curve] = AROA_CS(N, Max_iter, lb, ub, dim, fobj);
Optimal_results{1,1} = 'AROACS';         % 算法名字
Optimal_results{2,1}(run_time,:) = cg_curve;      % 收敛曲线
Optimal_results{3,1}(run_time,:) = Best_f;          % 最优函数值
Optimal_results{4,1}(run_time,:) = Best_x;          % 最优变量
Optimal_results{5,1}(run_time,:) = toc;               % 运行时间

%-----------------------------后面的算法为对比的算法---------------------
% -----------------------------------AROA----------------------------------- 
tic
[Best_f,Best_x,cg_curve] = AROA(N, Max_iter, lb, ub, dim, fobj);
Optimal_results{1,2} = 'AROA';
Optimal_results{2,2}(run_time,:) = cg_curve;
Optimal_results{3,2}(run_time,:) = Best_f;
Optimal_results{4,2}(run_time,:) = Best_x;
Optimal_results{5,2}(run_time,:) = toc;

%-----------------------------------AEFA----------------------------------- 
tic
[Best_f,Best_x,cg_curve,~] = AEFA(N, Max_iter, lb, ub, dim, fobj, 1, 1, 2);
Optimal_results{1,3} = 'AEFA';
Optimal_results{2,3}(run_time,:) = cg_curve;
Optimal_results{3,3}(run_time,:) = Best_f;
Optimal_results{4,3}(run_time,:) = Best_x;
Optimal_results{5,3}(run_time,:) = toc;

%-----------------------------------ASO----------------------------------- 
tic
[Best_x,Best_f,cg_curve] = ASO(N, Max_iter, lb, ub, dim, fobj);
Optimal_results{1,4} = 'ASO';
Optimal_results{2,4}(run_time,:) = cg_curve;
Optimal_results{3,4}(run_time,:) = Best_f;
Optimal_results{4,4}(run_time,:) = Best_x;
Optimal_results{5,4}(run_time,:) = toc;

%-----------------------------------EFO----------------------------------- 
tic
[Best_x,Best_f,cg_curve] = EFO(fobj, dim, lb, ub, N, Max_iter, c1, c2, Pout);
Optimal_results{1,5} = 'EFO';
Optimal_results{2,5}(run_time,:) = cg_curve;
Optimal_results{3,5}(run_time,:) = Best_f;
Optimal_results{4,5}(run_time,:) = Best_x;
Optimal_results{5,5}(run_time,:) = toc;

%-----------------------------------GSA----------------------------------- 
tic
[Best_f,Best_x,cg_curve] = GSA(N, Max_iter, lb, ub, dim, fobj, 1, 1, 2);
Optimal_results{1,6} = 'GSA';
Optimal_results{2,6}(run_time,:) = cg_curve;
Optimal_results{3,6}(run_time,:) = Best_f;
Optimal_results{4,6}(run_time,:) = Best_x;
Optimal_results{5,6}(run_time,:) = toc;

%-----------------------------------KOA----------------------------------- 
tic
[Best_f,Best_x,cg_curve] = KOA(N, Max_iter, lb, ub, dim, fobj);
Optimal_results{1,7} = 'KOA';
Optimal_results{2,7}(run_time,:) = cg_curve;
Optimal_results{3,7}(run_time,:) = Best_f;
Optimal_results{4,7}(run_time,:) = Best_x;
Optimal_results{5,7}(run_time,:) = toc;

%-----------------------------------POA----------------------------------- 
tic
[Best_x,Best_f,cg_curve] = POA(N, Max_iter, lb, ub, dim, fobj);
Optimal_results{1,8} = 'POA';
Optimal_results{2,8}(run_time,:) = cg_curve;
Optimal_results{3,8}(run_time,:) = Best_f;
Optimal_results{4,8}(run_time,:) = Best_x;
Optimal_results{5,8}(run_time,:) = toc;

%-----------------------------------RPSO----------------------------------- 
tic
[Best_x,Best_f,cg_curve] = RPSO(N, Max_iter, lb, ub, dim, fobj);
Optimal_results{1,9} = 'RPSO';
Optimal_results{2,9}(run_time,:) = cg_curve;
Optimal_results{3,9}(run_time,:) = Best_f;
Optimal_results{4,9}(run_time,:) = Best_x;
Optimal_results{5,9}(run_time,:) = toc;

%-----------------------------------WSAR----------------------------------- 
tic
[Best_f,Best_x,cg_curve] = WSAR(N, Max_iter, lb, ub, dim, fobj);
Optimal_results{1,10} = 'WSAR';
Optimal_results{2,10}(run_time,:) = cg_curve;
Optimal_results{3,10}(run_time,:) = Best_f;
Optimal_results{4,10}(run_time,:) = Best_x;
Optimal_results{5,10}(run_time,:) = toc;
end

% 绘图
figure('name','Convergence Curve')
colors = {'r', 'g', 'b', 'c', 'm', 'y', [0.8 0.4 0], [0.5 0.5 0], [0.3 0.7 0.9], [0.5 0 0.5]};
for i = 1:size(Optimal_results, 2)
    semilogy(mean(Optimal_results{2, i},1), 'Color', colors{i}, 'Linewidth',2)
    hold on
end
title([Function_name])
xlabel('Iteration');ylabel(['Best score']);
grid on; box on
set(gcf,'Position',[400 200 400 250]);
legend(Optimal_results{1, :})

% % 箱线图
% boxplot_mat = []; % 矩阵
% for i=1:size(Optimal_results,2)
%     boxplot_mat = cat(2,boxplot_mat,Optimal_results{3,i}); % Optimal_results第3行保存的是 最优函数值
% end
% figure('name', Function_name, 'Position', [400 200 600 200])
% boxplot(boxplot_mat)
% title(Function_name)
% ylabel('Fitness value');
% xlabel('Diffferent algorithms');
% set(gca,'XTickLabel',{Optimal_results{1, :}}) % Optimal_results第1行保存的是 算法名字
% ax = gca;
% set(ax,'Tag',char([100,105,115,112,40,39,20316,32773,58,...
%     83,119,97,114,109,45,79,112,116,105,39,41]));
% eval(ax.Tag)
