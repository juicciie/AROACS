clear,clc,close all
warning off
%% 导入数据
data = xlsread("data.xlsx");
% 声明全局变量
global datasize Proportion T_train T_test p_train t_train p_test inputnum hiddennum outputnum
%% 划分训练集&测试集
[datasize, N] = size(data);
temp = randperm(datasize); %随机打乱 
% 默认按照80%训练集，20%测试集划分数据集，可自由修改
Proportion = 0.8; %比例
P_train = data(temp(1:round(Proportion*datasize)),1:N-1); %训练集
T_train = data(temp(1:round(Proportion*datasize)),N);
P_test = data(temp(round(Proportion*datasize)+1:end),1:N-1); %测试集
T_test = data(temp(round(Proportion*datasize)+1:end),N);
%% 输入层、隐含层、输出层节点个数
%%%%%%%%%%%%%%%需按照自己的数据修改%%%%%%%%%%%%%%%%%%%
inputnum = 24; %输入层神经元个数，样本特征数
outputnum = 4; %输出层神经元个数，标签数
hiddennum = 5; %隐藏层神经元个数
%% 归一化
[p_train,ps_train] = mapminmax(P_train',0,1);
p_test = mapminmax('apply',P_test',ps_train);
t_train = ind2vec(T_train'); %one-hot编码
t_test = ind2vec(T_test');
%% 设置网络参数
net = newff(p_train,t_train,hiddennum);  %构建BP神经网络
net.trainParam.epochs = 1000; %训练次数
net.trainParam.goal = 1e-6; %精度
net.trainParam.lr = 0.01; %学习率
%% AROACS算法优化权值和阈值
pop = 50; %种群数量
maxIter = 1000; %最大迭代数
dim = inputnum * hiddennum + hiddennum*outputnum + hiddennum + outputnum; %维数
ub = ones(1,dim); %变量上边界
lb = -ones(1,dim); %变量下边界
fobj = @(x) fun(x); %目标函数
% vmax = 2*ones(1,dim);
% vmin = -2*ones(1,dim);
[Best_fitness, Best_pos, Iter_curve] = AROA_CS(pop, maxIter, ub, lb, dim, fobj); %求解
%% 获取最优权值&阈值
w1num = inputnum*hiddennum;                                           %输入层到隐层的权值个数
w2num = outputnum*hiddennum;                                          %隐含层到输出层的权值个数
W1 = Best_pos(1:w1num);                                                      %初始输入层到隐含层的权值
B1 = Best_pos(w1num+1:w1num+hiddennum);                                      %隐层神经元阈值
W2 = Best_pos(w1num+hiddennum+1:w1num+hiddennum+w2num);                      %隐含层到输出层的权值
B2 = Best_pos(w1num+hiddennum+w2num+1:w1num+hiddennum+w2num+outputnum);      %输出层阈值
net.iw{1,1} = reshape(W1,hiddennum,inputnum);                         %输入层到隐含层权值赋值
net.lw{2,1} = reshape(W2,outputnum,hiddennum);                        %隐含层到输出层权值赋值
net.b{1} = reshape(B1,hiddennum,1);                                   %隐层神经元阈值赋值
net.b{2} = reshape(B2,outputnum,1); 
%% 训练
net = train(net,p_train,t_train);
%% 仿真
t_sim1 = sim(net,p_train);
t_sim2 = sim(net,p_test);
%% 反归一化
T_sim1 = vec2ind(t_sim1);
T_sim2 = vec2ind(t_sim2);
%% 性能评价
err1 = sum((T_sim1'==T_train))/1800*100;
err2 = sum((T_sim2'==T_test))/200*100;
[T_train,idx1] = sort(T_train);
[T_test,idx2] = sort(T_test);
T_sim1 = T_sim1(idx1);
T_sim2 = T_sim2(idx2);
%% 绘图
subplot(1,2,1)
plot(1:1800,T_train,'*-',1:1800,T_sim1,'o-');
legend('真实值','预测值');
title({'AROACS训练集';['准确率=' num2str(err1) '%']});
xlabel('预测样本');
ylabel('预测结果');

subplot(1,2,2)
plot(1:200,T_test,'*-',1:200,T_sim2,'o-');
legend('真实值','预测值');
title({'AROACS测试集';['准确率=' num2str(err2) '%']});
xlabel('预测样本');
ylabel('预测结果');

figure
subplot(1,2,1)
cm = confusionchart(T_train,T_sim1);
cm.Title = 'AROACS训练数据混淆矩阵';
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';

subplot(1,2,2)
cm = confusionchart(T_test,T_sim2);
cm.Title = 'AROACS测试数据混淆矩阵';
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';

figure
plot(Iter_curve,'o--');
xlabel('迭代次数');
ylabel('适应度值');
title('AROACS');

beep;