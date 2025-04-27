% *************************************************************
% Evolutionary Field Optimization (EFO):
% Citation: Alagoz BB, Simsek OI, Ari D, Tepljakov A, Petlenkov E, Alimohammadi H.
% An Evolutionary Field Theorem: Evolutionary Field Optimization in Training 
% of Power-Weighted Multiplicative Neurons for Nitrogen Oxides-Sensitive 
% Electronic Nose Applications. Sensors. 2022; 22(10):3836. 
% https://doi.org/10.3390/s22103836
% *************************************************************
%% Configuration of optimization problem
% Objective function
fobj = @(C) SquareError(C); 
% Number of optimization parameters
dim = 30;
% Lower boundaries of optimization parameters 
MinValue=-100;
lb = MinValue*ones(1,dim);
% Upper boundaries of optimization parameters 
MaxValue=100;
ub = MaxValue*ones(1,dim);
%% Configuration of EFO
% Maximum Number of Iterations
MaxIterationNum = 100;
% Population Size
SearchAgentsNum = 20;         
% Crossover coefficient.
c1=0.26;
% Mutation coefficient.
c2=1.26;
%Absolute scattering size coefficient
Pout=0.31;
%Message flag.(0:stop messages. 1:Enables messages from optimization)1
MesFlag=1; 
%% Calling EFO
fprintf('EFO is starting.... \n');
%[Best_score,Best_pos,GWO_cg_curve]=GWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
[BestSol BestCost BestCostHist]=EFO(fobj,dim,lb,ub,SearchAgentsNum,MaxIterationNum,c1,c2,Pout,MesFlag);

%% Results from EFO
fprintf('Best Solution Vector: \n');
disp(BestSol);
fprintf('Best Cost: %d \n',BestCost);

figure(1)
semilogy(BestCostHist,'LineWidth', 2)
xlabel('Iteration')
ylabel('Square Error')
grid