function sol2=ParseSolution(sol1,model)
    % 从输入解中提取中间点坐标
    x = sol1(1:model.n);
    y = sol1(model.n+1:end);
    
    % 确保 x 和 y 的长度相同
    if length(x) ~= length(y)
        error('中间点的 x 和 y 坐标数量不一致');
    end
    
    % 构建完整路径点序列（包含起点和终点）
    xs=model.xs; % 起点
    ys=model.ys; 
    xt=model.xt; % 终点
    yt=model.yt;
    xobs=model.xobs; % 障碍物
    yobs=model.yobs;
    robs=model.robs;

    % 确保所有向量维度一致
    if ~isscalar(xs)
        xs = xs(1);
        ys = ys(1);
    end
    
    if ~isscalar(xt)
        xt = xt(1);
        yt = yt(1);
    end
    
    XS=[xs x xt]; % 完整x坐标序列
    YS=[ys y yt]; % 完整y坐标序列
    
    % 再次检查长度是否匹配
    if length(XS) ~= length(YS)
        error('XS和YS长度不匹配：XS长度=%d，YS长度=%d', length(XS), length(YS));
    end
    
    k=numel(XS); % 路径点数量
    TS=linspace(0,1,k); % 生成均匀时间序列
    
    tt=linspace(0,1,100); % 更密的时间序列
    xx=spline(TS,XS,tt); % x方向样条插值
    yy=spline(TS,YS,tt); % y方向样条插值
    
    % 计算路径长度
    dx=diff(xx); % x方向差分
    dy=diff(yy); % y方向差分
    
    L=sum(sqrt(dx.^2+dy.^2)); % 总路径长度
    
    nobs = numel(xobs); % 障碍物的数量
    Violation = 0;
    for k=1:nobs
        d=sqrt((xx-xobs(k)).^2+(yy-yobs(k)).^2); % 计算到障碍物的距离
        v=max(1-d/robs(k),0); % 计算违反程度
        Violation=Violation+mean(v);
    end
    
    sol2.TS=TS;
    sol2.XS=XS;
    sol2.YS=YS;
    sol2.tt=tt;
    sol2.xx=xx;
    sol2.yy=yy;
    sol2.dx=dx;
    sol2.dy=dy;
    sol2.L=L;
    sol2.Violation=Violation;
    sol2.IsFeasible=(Violation==0);
end