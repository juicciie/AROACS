function [fbest, xbest, Convergence_curve, St] = AROA_CS(N, maxEvals, lb, ub, dim, fobj)
    
    % Algorithm parameters definition
    c = 0.98;      % 增大到0.98,增强收敛能力
    fr1 = 0.1;     % 减小到0.1,用于局部精细搜索
    fr2 = 0.4;     % 减小到0.4,减少随机性
    p1 = 0.3;      % 增大到0.3,增加向最优解学习的概率  
    p2 = 0.7;      % 适当减小,平衡全局与局部搜索
    Ef = 0.3;      % 减小扰动概率
    tr1 = 0.95;    % 增大,提高局部搜索精度
    tr2 = 0.9;     % 增大,增强核心区域搜索
    tr3 = 0.85;    % 略微降低,保持适度探索
    % Algorithm parameters definition

    % New parameters for anti-attraction
    alpha = 0.15;   % 较小的反吸引强度
    beta = 0.3;    % 较小的适应度差异阈值
    
    tmax = maxEvals;
    maxEvals = 2 * N *maxEvals;
    evalCounter = 0;

    Convergence_curve = zeros(1,tmax);
    Xmin = repmat(ones(1,dim).*lb,N,1);
    Xmax = repmat(ones(1,dim).*ub,N,1);

    % random initialization - Eq (3)
    % 改进策略1:Hénon-Sine hyperchaotic initialization
    X = henon_sine_initialization(N, dim, lb, ub);
    [X, F, evalCounter] = evaluate_population(X, fobj, ub, lb, evalCounter, maxEvals);
    [fbest, ibest] = min(F);
    xbest = X(ibest,:);
    % random initialization - Eq (3)

    X_memory = X;
    F_memory = F;

    % 添加历史最优解记录
    history_best = zeros(tmax, dim);
    history_best(1,:) = xbest;

    % Main loop
    for t=1:tmax	
        D = squareform(pdist(X, 'squaredeuclidean'));  % Eq (4) 
        m = tanh(t, tmax, [-2, 7]);   % Eq (11) 

        % 自适应控制参数
        c_adaptive = c * (1 - 0.5*t/tmax);

        for i=1:N
           Dimax = max(D(i,:));
           k = floor((1-t/tmax)*N)+1;  % Eq (9)
           [~, neighbors] = sort(D(i,:));
           
           % Attraction-Repulsion operator % Eq (6)
           % 改进策略2:反吸引策略
           delta_ni = zeros(1,dim);
           for j = neighbors(1:k)
                I = 1 - (D(i,j) / Dimax);  % Eq (7)
                % 计算适应度差异
                fitness_diff = abs(F(j)-F(i)) / max(abs(F));

                if fitness_diff > beta
                    % 反吸引策略
                    s = -sign(F(j)-F(i)) * alpha * (1 + fitness_diff);
                else
                    % 原始吸引策略
                    s = sign(F(j)-F(i));  % Eq (8)
                end

                % 添加扰动项
                disturb = 0.1 * rand() * (1-t/tmax);        
                delta_ni = delta_ni + c_adaptive*(X_memory(i,:) - X_memory(j,:)) * I * s * (1 + disturb);
           end
           ni = delta_ni/N;
           % Attraction-Repulsion operator % Eq (6)

           % Attraction to best solusion Eq (10)
           % Enhanced attraction to best solution
           if rand < p1
               bi = m*c_adaptive.*(rand(1,dim).*xbest - X_memory(i,:));
           else
               % 添加历史信息指导
               if t > 1
                   prev_best = history_best(t-1,:);
                   bi = m*c_adaptive.*(0.7*xbest + 0.3*prev_best - X_memory(i,:));
               else
                   bi = m*c_adaptive.*(xbest - X_memory(i,:));
               end
           end
           % Attraction to best solusion Eq (10)

           % Local search operators Eq (15)
           % 改进策略3:核心空间定位
           if rand < p2
               if rand > 0.5 * t / tmax + 0.25
                   % 核心空间定位策略
                   core_radius = (ub - lb) * (1 - t / tmax)^1.5; % 核心区域半径
                   core_center = xbest; % 以当前最优解为核心中心
                   
                   % 计算当前解到核心空间的距离
                   dist_to_core = sqrt(sum((X(i,:)-core_center).^2));

                   % 在核心区域内
                   if dist_to_core <= norm(core_radius)
                       % 精细搜索
                       u1 = rand(1, dim) > tr1;
                       ri = u1.*random('Normal', zeros(1,dim), fr1*(1-t/tmax)*(ub-lb)/2); % 减小搜索范围
                   else
                        u1 = rand(1, dim) > tr1;
                        ri = u1.*random('Normal', zeros(1,dim), fr1*(1-t/tmax)*(ub-lb));
                   end

                   % 添加定向搜索分量
                   direction_to_core = (core_center - X(i,:))/dist_to_core;
                   ri = ri + 0.3 * (1 - t / tmax) * direction_to_core; % 0.3为定向强度系数

               else
                   % 基于核心区域的自适应搜索
                   u2 = rand(1,dim) > tr2;
                   w = index_roulette_wheel_selection(F, k);
                   Xw = X_memory(w,:);

                   % 计算搜索步长
                   dist_factor = min(1, sqrt(sum((X(i,:) - xbest).^2))/(norm(ub-lb)));
                   search_scale = fr2*(1-t/tmax)*dist_factor;

                   if rand < 0.5
                       ri = search_scale*u2.*sin(2*pi*rand(1,dim)).*abs(rand(1,dim).*Xw-X_memory(i,:));
                   else
                       ri = search_scale*u2.*cos(2*pi*rand(1,dim)).*abs(rand(1,dim).*Xw-X_memory(i,:));
                   end

                   % 添加核心吸引
                   core_attraction = 0.15*(1-t/tmax)*(xbest - X(i,:));
                   ri = ri + core_attraction;
               end
           else
               % 全局探索
               u3 = rand(1,dim) > tr3;
               ri = u3.*(2*rand(1,dim)-ones(1,dim)) .* (ub-lb);

               % 添加轻微的核心引导
               if rand < 0.3 % 30%概率进行核心引导
                   ri = ri + 0.05 * (xbest - X(i,:));
               end
           end
           % Local search operators Eq (15)

           X(i,:) = X(i,:) + ni + bi + ri;  % Eq(16)
        end
        
        [X, F, evalCounter] = evaluate_population(X, fobj, ub, lb, evalCounter, maxEvals);
        [fbest_candidate, ibest_candidate] = min(F);

        if fbest_candidate < fbest
            fbest = fbest_candidate;
            xbest = X(ibest_candidate, :);
        end

        % 更新历史最优解
        history_best(t,:) = xbest;

	    [X, F] = memory_operator(X, F, X_memory, F_memory);  % Eq (18)
        X_memory = X;
        F_memory = F; 
			
        % Eq (17)			 
        CF=(1-t/tmax)^3;			 
        if rand < Ef
            u4 = rand(N,dim) < Ef;                                                                                              
            X = X + CF*(u4.*(rand(N,dim).*(Xmax-Xmin) + Xmin));
        else
            r7 = rand();		
            X = X + (CF*(1-r7) + r7)*(X(randperm(N),:) - X(randperm(N),:));
        end
		% Eq (17)

        [X, F, evalCounter] = evaluate_population(X, fobj, ub, lb, evalCounter, maxEvals);
        [fbest_candidate, ibest_candidate] = min(F);

        if fbest_candidate < fbest
            fbest = fbest_candidate;
            xbest = X(ibest_candidate, :);
        end
	 
        [X, F] = memory_operator(X, F, X_memory, F_memory);  % Eq (18)
        X_memory = X;
        F_memory = F;

  	    Convergence_curve(t) = fbest;
    end

    for i = 1 : N    
        [fit(i), st(i)] = fobj(X(i, :)) ;              
    end
    [ ~, bestI ] = min( fit ); 
    St = st(bestI);
		
end

function [X] = henon_sine_initialization(N, dim, lb, ub)
    % Two-dimensional Hénon-Sine hyperchaotic map initialization
    
    % Initialize population matrix
    X = zeros(N, dim);
    
    % Handle scalar bounds
    if length(lb) == 1
        lb = lb * ones(1,dim);
    end
    if length(ub) == 1 
        ub = ub * ones(1,dim);
    end
    
    % Hénon-Sine map parameters
    a = 1.4;  % control parameter
    b = 0.3;  % control parameter
    
    for i = 1:N
        % Random initial values in [0,1]
        x = rand();
        y = rand();
        
        % Generate sequence for each dimension
        for j = 1:dim
            % Hénon-Sine map iterations
            x_new = 1 - a*(sin(pi*x))^2 + y;
            y_new = b*x;
            
            % Update values
            x = x_new;
            y = y_new;
            
            % Map to search space and store
            % Using x value for initialization as it has better ergodicity
            X(i,j) = mod(abs(x), 1) * (ub(j)-lb(j)) + lb(j);
        end
    end
    
    % Add small perturbation to enhance diversity
    X = X + 0.005 * randn(size(X)) .* (ub-lb);
    
    % Ensure bounds
    X = max(min(X, ub), lb);
end

function [X, F, evalCounter] = evaluate_population(X, fobj, ub, lb, evalCounter, maxEvals)
    N = size(X,1);
    F = Inf(N,1);
    X = max(lb, min(ub, X)); % Check space bounds
    
    for i=1:N
        if evalCounter >= maxEvals
            break
        end
        F(i) = fobj(X(i,:));
        evalCounter = evalCounter + 1;
    end
end

function [X, F] = memory_operator(X, F, X_memory, F_memory)
    dim = size(X, 2);
    Inx = F_memory < F;
    Indx = repmat(Inx,1,dim);
    X = Indx.*X_memory + ~Indx.*X;
    F = Inx.*F_memory + ~Inx.*F;
end

function [y] = tanh(t, tmax, range)
    z = 2*(t/tmax*(range(2)-range(1)) + range(1));
    y = 0.5*((exp(z)-1)/(exp(z)+1) + 1);
end


function [selected_index] = index_roulette_wheel_selection(F, k)
    fitness = F(1:k);
    weights = max(fitness) - fitness;
    weights = cumsum(weights/sum(weights));
    
    selected_index = roulette_wheel_selection(weights);
end

function [selected_index] = roulette_wheel_selection(weights)
    r = rand();
    selected_index = 1;
    for index=size(weights,1)
        if r <= weights(index)
            selected_index = index;
            break;
        end
    end
end

% 在全局搜索阶段添加最有希望抽样策略
function [X, F] = promising_sampling_strategy(X, F, lb, ub, dim, fobj, evalCounter, maxEvals)
    % 计算种群的适应度分布
    fitness_prob = (max(F) - F) ./ (max(F) - min(F));
    fitness_prob = fitness_prob ./ sum(fitness_prob);
    
    % 选择最有希望的个体
    num_promising = floor(size(X, 1) * 0.2);  % 选择20%最有希望的个体
    [~, promising_indices] = maxk(fitness_prob, num_promising);
    
    % 对最有希望的个体进行扩展
    promising_solutions = X(promising_indices, :);
    
    % 使用高斯变异生成新解
    for i = 1:num_promising
        % 高斯扰动
        mutation_strength = 0.1 * (1 - evalCounter / maxEvals);
        new_solution = promising_solutions(i, :) + ...
            mutation_strength * randn(1, dim) .* (ub - lb);
        
        % 边界处理
        new_solution = max(min(new_solution, ub), lb);
        
        % 评估新解
        new_fitness= fobj(new_solution);
        
        % 如果新解优于最差解，则替换
        [worst_fitness, worst_idx] = max(F);
        if new_fitness < worst_fitness
            X(worst_idx, :) = new_solution;
            F(worst_idx) = new_fitness;
        end
    end
    
    % 重新评估种群
    [X, F, evalCounter] = evaluate_population(X, fobj, ub, lb, evalCounter, maxEvals);
end




