function [fbest, xbest, Convergence_curve] = AROA_CS(N, maxEvals, lb, ub, dim, fobj)
    
    % Algorithm parameters definition
    c = 0.85;      % 降低以减少过度收敛 (原0.95)
    fr1 = 0.12;    % 减小局部搜索范围 (原0.15)
    fr2 = 0.45;    % 降低以增加精细搜索 (原0.6)
    p1 = 0.25;     % 略微提高最优解吸引概率 (原0.2)
    p2 = 0.75;     % 略微降低局部搜索概率 (原0.8)
    Ef = 0.35;     % 降低扰动强度 (原0.4)
    tr1 = 0.92;    % 提高精细搜索触发概率 (原0.9)
    tr2 = 0.88;    % 提高中等范围搜索概率 (原0.85)
    tr3 = 0.85;    % 降低全局搜索概率 (原0.9)
    
    % 反吸引参数
    alpha = 0.25;   % 降低反吸引强度 (原0.3)
    beta = 0.4;     % 降低适应度差异阈值 (原0.5)

    tmax = maxEvals;
    maxEvals = 2 * N *maxEvals;
    evalCounter = 0;

    Convergence_curve = zeros(1,tmax);
    Xmin = repmat(ones(1,dim).*lb,N,1);
    Xmax = repmat(ones(1,dim).*ub,N,1);


    % random initialization - Eq (3)
    % Hénon-Sine hyperchaotic initialization
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
           % 改进策略1:反吸引策略
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
           % 改进策略2:核心空间定位
           if rand < p2
               if rand > 0.5 * t / tmax + 0.25
                   % 核心空间定位策略
                   core_radius = (ub - lb) * (1 - t / tmax); % 核心区域半径 
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
                   core_attraction = 0.1*(1-t/tmax)*(xbest - X(i,:));
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
end

function [X] = henon_sine_initialization(N, dim, lb, ub)
    % Two-dimensional Hénon-Sine hyperchaotic map initialization
    
    % Initialize population matrix
    X = zeros(N, dim);
    
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
            X(i,j) = mod(abs(x), 1) * (ub-lb) + lb;
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
