clear
clc
close all
addpath input_data

Max_iter = 2000; 
N = 50;
runs = 5;
Function_name = 'F18';
path = '/Users/juicciie/Paper/complexity evaluation/100D/';
elapsed_time = {};
[lb, ub, dim, fobj] = CEC2017_100(Function_name);
for i = 1:runs
    tic;
    [Best_f,Best_x,cg_curve] = AROA(N, Max_iter, lb, ub, dim, fobj);
    elapsed_time{i} = toc;
end

save([path, 'AROA'], 'elapsed_time');