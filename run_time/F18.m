% 计算CEC测试集F18函数的运行时间
clear all;
clc;

% 设置问题维度
dim = 100; % 可以根据需要修改维度

% 设置搜索空间边界（根据CEC测试集的定义）
lb = -100 * ones(1, dim);
ub = 100 * ones(1, dim);

% 随机生成一个解进行测试
x = initialization(50, dim, ub, lb);

% 记录开始时间
tic;

% 调用F18函数进行评估
% 注意：根据不同版本的CEC测试集，函数名和参数可能略有不同
fitness = cec17_func(x, 18); % 第二个参数18表示F18函数

% 记录结束时间并计算总用时
elapsed_time = toc;

% 显示结果
fprintf('F18函数评估值: %e\n', fitness);
fprintf('运行时间: %f 秒\n', elapsed_time);