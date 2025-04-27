% ���������Ż��㷨 
function [Best_pos,Best_fitness,Iter_curve, History_pos, History_best]=GSA(pop, maxIter, lb, ub, dim,fobj)
%input
%pop ��Ⱥ����
%dim ����ά��
%ub �����ϱ߽�
%lb �����±߽�
%fobj ��Ӧ�Ⱥ���
%maxIter ����������
%output
%Best_pos ����λ��
%Best_fitness ������Ӧ��ֵ
%Iter_curve ÿ��������Ӧ��ֵ
%History_pos ÿ����Ⱥλ��
%History_best ÿ�����Ÿ���λ��
%% ��ʼ����Ⱥλ��
X = initialization(pop, ub, lb, dim);
%% ������Ӧ��ֵ
for i = 1:pop
    fitness(i) = fobj(X(i,:));
end
% ����λ��&������Ӧ��ֵ
[SortFitness, indexSort] = sort(fitness);
gBest = X(indexSort(1),:);
gBestFitness = SortFitness(1);
M = zeros(pop, 1); %��������
V = zeros(pop,dim); %�ٶȾ���
%% ����
for t = 1:maxIter
    %��������
    M = massCalculation(fitness);
    %������������
    G = Gconstant(t, maxIter);
    %������ٶ�
    a = Acceleration(M,X,G,t,maxIter);
    %λ�ø���
    [X,V] = move(X,a,V);
    %�߽���
    Flag4ub=X(i,:)>ub;
    Flag4lb=X(i,:)<lb;
    X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    
    %������Ӧ��ֵ
    for i = 1:pop
        fitness(i) = fobj(X(i,:));
        if fitness(i) < gBestFitness
            gBestFitness = fitness(i);
            gBest = X(i,:);
        end
    end
    History_pos{t} = X;
    History_best{t} = gBest;
    Iter_curve(t) = gBestFitness;
end
Best_pos = gBest;
Best_fitness = gBestFitness;
end
%% ��ʼ������
function X=initialization(SearchAgents_no,ub,lb,dim)

Boundary_no= size(ub,2); 
if Boundary_no==1
    X=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end
end
%% ��������
function M = massCalculation(fitness)
    bestF = min(fitness);
    worstF = max(fitness);
    M = (fitness - worstF) ./ (bestF - worstF);
    M = M ./ sum(M);
end
%% ������������
function G = Gconstant(iter, max_it)
    alfa = 20;
    G0 = 100;
    G = G0*exp(-alfa*iter/max_it);
end
%% ������ٶ�
function a = Acceleration(M, X, G, iter, max_it)
    [N, dim] = size(X);
    final_per = 2;
    kbest = final_per + (1 - iter/max_it) * (100-final_per);
    kbest = round(N*kbest/100);
    [Ms, ds] = sort(M, 'descend');
    for i = 1:N
        F(i,:) = zeros(1,dim);
        for ii = 1:kbest
            j = ds(ii);
            if j ~= i
                R = norm(X(i,:) - X(j,:), 2);
                for k = 1:dim
                    F(i,k) = F(i,k) + rand*M(j)*((X(j,k) - X(i,k))/(R+eps));
                end
            end
        end
    end
    a = F.*G;
end
%% λ�ø���
function [X,V] = move(X,a,V)
    [N,dim] = size(X);
    V = rand(N,dim).*V+a;
    X = X + V;
end