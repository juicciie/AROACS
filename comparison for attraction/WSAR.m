function [fbest, xbest, Convergence_curve] = WSAR(N, maxEvals, lb, ub, dim, fobj)
% WSAR: Weighted Superposition Attraction/Repulsion Algorithm 
% Prof. Dr. Adil Baykasoğlu, DEU IE, İzmir, Türkiye
% 
% 输入参数:
%   N          - 种群规模
%   maxEvals   - 最大评估次数
%   lb         - 下界向量
%   ub         - 上界向量
%   dim        - 问题维度
%   fobj       - 目标函数句柄
% 
% 输出参数:
%   fbest            - 找到的最优解的目标函数值
%   xbest            - 找到的最优解
%   Convergence_curve - 收敛曲线（每次评估的最优值）

    % 算法参数
    tao = -0.8;  % 权重参数

    % 确保lb和ub是向量
    if length(lb) == 1
        lb = lb * ones(1, dim);
    end
    
    if length(ub) == 1
        ub = ub * ones(1, dim);
    end
    
    % 初始化种群
    agent = zeros(N, dim);
    for i = 1:N
        for j = 1:dim
            agent(i, j) = lb(j) + rand() * (ub(j) - lb(j));
        end
    end
    
    % 计算初始目标函数值
    agent_OBJ = zeros(N, 1);
    
    % 初始化收敛曲线 - 现在大小为maxEvals
    Convergence_curve = zeros(1, maxEvals);
    evaluations = 0;
    
    % 记录最优解
    fbest = Inf;
    xbest = [];
    
    % 计算初始种群的目标函数值
    for i = 1:N
        evaluations = evaluations + 1;
        agent_OBJ(i) = fobj(agent(i, :));
        
        % 更新全局最优
        if agent_OBJ(i) < fbest
            fbest = agent_OBJ(i);
            xbest = agent(i, :);
        end
        
        % 记录每次评估后的最优值
        Convergence_curve(evaluations) = fbest;
    end
    
    % 主循环
    while evaluations < maxEvals
        % 根据目标函数值排序
        [~, sortIdx] = sort(agent_OBJ);
        index = zeros(1, N);
        for i = 1:N
            index(sortIdx(i)) = i;
        end
        reverse_index = N - index + 1;
        
        % 计算权重
        BESTweight = index.^tao;
        WORSTweight = reverse_index.^tao;
        
        % 随机向量
        randomvector = rand(1, dim);
        
        % 确定最佳超位置
        BESTsuperposition = zeros(1, dim);
        for i = 1:dim
            pos_candidate1 = [];
            p1 = [];
            for j = 1:N
                if BESTweight(j) >= randomvector(i)
                    pos_candidate1 = [pos_candidate1; agent(j, i)];
                    p1 = [p1; BESTweight(j)];
                end
            end
            
            if ~isempty(pos_candidate1)
                % 使用权重抽样
                cumProbs = cumsum(p1) / sum(p1);
                r = rand();
                idx = find(cumProbs >= r, 1, 'first');
                BESTsuperposition(i) = pos_candidate1(idx);
            else
                % 如果没有满足条件的候选者
                BESTsuperposition(i) = agent(randi(N), i);
            end
        end
        
        BEST_superposition = BESTsuperposition;
        evaluations = evaluations + 1;
        BESTsuperposition_OBJ = fobj(BEST_superposition);
        
        % 更新全局最优
        if BESTsuperposition_OBJ < fbest
            fbest = BESTsuperposition_OBJ;
            xbest = BEST_superposition;
        end
        
        % 记录每次评估后的最优值
        Convergence_curve(evaluations) = fbest;
        
        if evaluations >= maxEvals
            break;
        end
        
        % 确定最差超位置
        WORSTsuperposition = zeros(1, dim);
        for i = 1:dim
            pos_candidate2 = [];
            p2 = [];
            for j = 1:N
                if WORSTweight(j) >= randomvector(i)
                    pos_candidate2 = [pos_candidate2; agent(j, i)];
                    p2 = [p2; WORSTweight(j)];
                end
            end
            
            if ~isempty(pos_candidate2)
                % 使用权重抽样
                cumProbs = cumsum(p2) / sum(p2);
                r = rand();
                idx = find(cumProbs >= r, 1, 'first');
                WORSTsuperposition(i) = pos_candidate2(idx);
            else
                % 如果没有满足条件的候选者
                WORSTsuperposition(i) = agent(randi(N), i);
            end
        end
        
        WORST_superposition = WORSTsuperposition;
        evaluations = evaluations + 1;
        WORSTsuperposition_OBJ = fobj(WORST_superposition);
        
        % 更新全局最优
        if WORSTsuperposition_OBJ < fbest
            fbest = WORSTsuperposition_OBJ;
            xbest = WORST_superposition;
        end
        
        % 记录每次评估后的最优值
        Convergence_curve(evaluations) = fbest;
        
        if evaluations >= maxEvals
            break;
        end
        
        % 如果需要，交换超位置
        if WORSTsuperposition_OBJ < BESTsuperposition_OBJ
            dummyx = BEST_superposition;
            BEST_superposition = WORST_superposition;
            WORST_superposition = dummyx;
            
            dummyOBJ = BESTsuperposition_OBJ;
            BESTsuperposition_OBJ = WORSTsuperposition_OBJ;
            WORSTsuperposition_OBJ = dummyOBJ;
        end
        
        % 生成新解
        nagent = zeros(N, dim);
        
        for j1 = 1:N
            if evaluations >= maxEvals
                break;
            end
            
            nagent(j1, :) = agent(j1, :);
            
            if agent_OBJ(j1) > BESTsuperposition_OBJ
                nagent(j1, :) = agent(j1, :) + (rand() * (BEST_superposition - abs(agent(j1, :))) + rand() * (abs(agent(j1, :)) - WORST_superposition));
            end
            
            if agent_OBJ(j1) <= BESTsuperposition_OBJ
                nagent(j1, :) = random_move(nagent(j1, :), dim);
            end
            
            % 边界处理
            for d = 1:dim
                if nagent(j1, d) < lb(d)
                    nagent(j1, d) = lb(d);
                end
                if nagent(j1, d) > ub(d)
                    nagent(j1, d) = ub(d);
                end
            end
            
            % 计算新解的目标函数值
            evaluations = evaluations + 1;
            new_obj = fobj(nagent(j1, :));
            
            % 更新全局最优
            if new_obj < fbest
                fbest = new_obj;
                xbest = nagent(j1, :);
            end
            
            % 记录每次评估后的最优值
            Convergence_curve(evaluations) = fbest;
            
            % 择优选择
            if new_obj <= agent_OBJ(j1)
                agent(j1, :) = nagent(j1, :);
                agent_OBJ(j1) = new_obj;
            end
        end
    end
    
    % 裁剪收敛曲线到实际评估次数
    Convergence_curve = Convergence_curve(1:evaluations);
end

% 随机移动函数
function x = random_move(x, dim)
    step_size = rand();
    for i = 1:dim
        x(i) = x(i) + (rand() * 2 - 1) * step_size;
    end
end