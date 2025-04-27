% RPSO - Repulsive Particle Swarm Optimization with Visualization

function [bestSol, bestFitness, convergence_curve] = RPSO(nPop, maxIter, lb, ub, dim, fobj)
    % 参数初始化
    w = 0.7;              % 惯性权重
    c1 = 1.5;             % 个体学习因子
    c2 = 1.5;             % 社会学习因子
    repulsion_coeff = 0.1; % 排斥系数
    bestFitness = inf;     % 全局最优适应度
    bestSol = zeros(1, dim); % 全局最优解

    % 初始化粒子的位置和速度
    particles = lb + (ub-lb) .* rand(nPop, dim); % 粒子的位置
    velocities = zeros(nPop, dim);              % 粒子的速度
    pBest = particles;                          % 个体最优位置
    pBestFitness = inf(nPop, 1);                % 个体最优适应度
    convergence_curve = zeros(maxIter, 1);      % 收敛曲线数据

    % 主循环
    for iter = 1:maxIter
        % 计算适应度
        fitness = zeros(nPop, 1);
        for i = 1:nPop
            fitness(i) = fobj(particles(i, :));
            % 更新个体最优
        end

        % 更新全局最优
        [minFitness, minIndex] = min(fitness);
        if minFitness < bestFitness
            bestFitness = minFitness;
            bestSol = particles(minIndex, :);
        end

        % 记录收敛曲线
        convergence_curve(iter) = bestFitness;

        % 更新速度和位置
        for i = 1:nPop
            r1 = rand(1, dim);
            r2 = rand(1, dim);

            % 排斥项的计算
            repulsionForce = zeros(1, dim);


            % 速度更新公式 (包括排斥项)
            velocities(i, :) = w * velocities(i, :) ...
                             + c1 * r1 .* (pBest(i, :) - particles(i, :)) ...
                             + c2 * r2 .* (bestSol - particles(i, :)) ...
                             + repulsionForce;

            % 位置更新
            particles(i, :) = particles(i, :) + velocities(i, :);

            % 边界控制
        end
    end
end
