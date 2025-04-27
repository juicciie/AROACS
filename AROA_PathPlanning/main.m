clear,clc,close all
rng(1)

%% 地图参数
model = CreateModel(); %构建地图模型
model.n = 3; %中间结点

%障碍物参数
model.xobs=[1.5 4.0 1.2 5 5];
model.yobs=[4.5 3.0 1.5 1 5];
model.robs=[1 1 1 1 1 1];

%范围限制
model.xmin = -10; model.xmax = 10;
model.ymin = -10; model.ymax = 10;

%起点
model.xs = 0;
model.ys = 0;
%终点
model.xt = 6;
model.yt = 6;

%% 算法参数
N = 50; %种群数量
Max_iter = 1000; %迭代次数
run_times = 30;
dim = model.n*2; %维数

% EFO default configurations
c1=0.26;
c2=1.26;
Pout=0.31;

lb = [model.xmin*ones(1,model.n) model.ymin*ones(1,model.n)];
ub = [model.xmax*ones(1,model.n) model.ymax*ones(1,model.n)];
fobj = @(x) MyCost(x,model);
Optimal_results={}; 
%% 求解
for run_time = 1:run_times
% -----------------------------------你的算法放在首位--------------------------------
tic
[Best_f,Best_x,cg_curve,ST1] = AROA_CS(N, Max_iter, lb, ub, dim, fobj);
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
[Best_f,Best_x,cg_curve] = GSA(N, Max_iter, lb, ub, dim, fobj);
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
%% 绘图
PlotSolution(ST1,model);

figure('name','Convergence Curve')
colors = {'r', 'g', 'b', 'c', 'm', 'y', [0.8 0.4 0], [0.5 0.5 0], [0.3 0.7 0.9], [0.5 0 0.5]};
for i = 1:size(Optimal_results, 2)
    semilogy(mean(Optimal_results{2, i},1), 'Color', colors{i}, 'Linewidth',2)
    hold on
end
title(['Path Planning'])
xlabel('Iteration');ylabel(['Best score']);
grid on; box on
set(gcf,'Position',[400 200 400 250]);
legend('AROACS', 'AROA', 'AEFA', 'ASO', 'EFO', 'GSA', 'KOA', 'POA','RPSO', 'WSAR')