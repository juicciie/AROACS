% 代价函数，用于评估路径的质量
% 代价函数考虑两个因素：路径长度(L)和障碍物违反程度(Violation)
function [z, sol]=MyCost(sol1,model)
        
    sol=ParseSolution(sol1,model); % 首先解析解决方案
    
    beta=100; % 惩罚因子，用于加重对障碍物冲突的惩罚
    % 如果路径没有与障碍物冲突(Violation=0)，总代价就等于路径长度
    z=sol.L*(1+beta*sol.Violation); % 计算总代价
    % 如果有冲突，代价会显著增加，引导优化算法寻找无冲突路径

end