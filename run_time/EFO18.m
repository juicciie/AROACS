clear
clc
close all
addpath input_data

Max_iter = 4000; 
N = 50;
runs = 5;
Function_name = 'F18';
c1 = 0.26;
c2 = 1.26;
Pout = 0.31;
path = '/Users/juicciie/Paper/complexity evaluation/100D/';
elapsed_time = {};
[lb, ub, dim, fobj] = CEC2017_100(Function_name);
for i = 1:runs
    tic;
    [Best_x,Best_f,cg_curve] = EFO(fobj, dim, lb, ub, N, Max_iter, c1, c2, Pout);
    elapsed_time{i} = toc;
end

save([path, 'EFO'], 'elapsed_time');