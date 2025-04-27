function fitness = fun(x)
global datasize Proportion T_train T_test p_train t_train p_test inputnum hiddennum outputnum

%% 构建BP神经网络
net=newff(p_train,t_train,hiddennum);  %构建BP神经网络
net.trainParam.epochs=1000; %训练次数
net.trainParam.goal=1e-6; %精度
net.trainParam.lr=0.01; %学习率
net.trainParam.showwindow=false; %关闭图形框
%% BP神经网络初始权值和阈值
w1num=inputnum*hiddennum;                                           %输入层到隐层的权值个数
w2num=outputnum*hiddennum;                                          %隐含层到输出层的权值个数
W1=x(1:w1num);                                                      %初始输入层到隐含层的权值
B1=x(w1num+1:w1num+hiddennum);                                      %隐层神经元阈值
W2=x(w1num+hiddennum+1:w1num+hiddennum+w2num);                      %隐含层到输出层的权值
B2=x(w1num+hiddennum+w2num+1:w1num+hiddennum+w2num+outputnum);      %输出层阈值
net.iw{1,1}=reshape(W1,hiddennum,inputnum);                         %为神经网络的输入层到隐含层权值赋值
net.lw{2,1}=reshape(W2,outputnum,hiddennum);                        %为神经网络的隐含层到输出层权值赋值
net.b{1}=reshape(B1,hiddennum,1);                                   %为神经网络的隐层神经元阈值赋值
net.b{2}=reshape(B2,outputnum,1);                                   %为神经网络的输出层阈值赋值
%% 训练
net = train(net,p_train,t_train);
%% 仿真
t_sim1 = sim(net,p_train);
t_sim2 = sim(net,p_test);
%% 反归一化
T_sim1 = vec2ind(t_sim1);
T_sim2 = vec2ind(t_sim2);
%% 性能评价
err1 = sum((T_sim1'~=T_train))/(datasize*Proportion);
err2 = sum((T_sim2'~=T_test))/(datasize*(1-Proportion));

% 训练集和测试集误差
fitness = err1 + err2;
end