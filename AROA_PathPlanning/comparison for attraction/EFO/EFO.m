% *****************************************
% Evolutionary Field Optimization (EFO):
% Citation: Alagoz BB, Simsek OI, Ari D, Tepljakov A, Petlenkov E, Alimohammadi H.
% An Evolutionary Field Theorem: Evolutionary Field Optimization in Training 
% of Power-Weighted Multiplicative Neurons for Nitrogen Oxides-Sensitive 
% Electronic Nose Applications. Sensors. 2022; 22(10):3836. 
% https://doi.org/10.3390/s22103836
%*******************************************************************
% Input Paremeters:
% fobj: Objective function (Defined in the form of 
% fobj = @(X) ObjectiveFunctionName(X))
% dim: Dimension of search space (Number of component in the property code X 
% that search agents are represented in the solution environment of the problem. 
% Lb: Vector that defines lower boundaries of the property code X in the search space.
% Ub: Vector that defines upper boundaries of the property code X in the search space.
% NumOfAgents: Number of search agents that are used in the solution environment.
% (Population size of agents)
% c1: Crossover coefficient. 
% c2: Mutation coefficient. 
% Pout: Absolute scattering size coefficent
% MesFlag: Message flag.(0 for no-message during optimization and 1 for message during optimization) 
%****************************************************************
% Output Parameters:
% BestSol: The optimal solution of EFO algorithim (The best the property code X)
% BestCost: The cost value of the optimal solution BestSol (The objective
% function value F(X)
% BestCostHist: Contains the best cost values during EFO optimization
%****************************************************************

function [BestSol, BestCost, BestCostHist]=EFO(fobj,dim,Lb,Ub,NumOfAgents,NumOfIteration,c1,c2,Pout)

for agent=1:NumOfAgents
    if dim == length(Lb) && dim == length(Ub)
    for var=1:dim
        L(agent,var)=Lb(var);
        U(agent,var)=Ub(var);
    end
    elseif length(Lb)==1 || length(Ub)==1
    for var=1:dim
        L(agent,var)=Lb;
        U(agent,var)=Ub;
    end
    else
      fprintf('Error : size of search range vector (lower and upper boundary) should be equal to number of the variable! \n'); 
    end
end
% Assigns initial value of agent property vector (X) randomly 
% between upper boundry (Ub) and lower boundray (Lb)
for agent=1:NumOfAgents
    for var=1:dim
        X(agent,var)=L(agent,var)+(U(agent,var)-L(agent,var))*rand;
    end
end
Xn=X;
% Calculates field values of initial agents (F(X)).
for agent=1:size(X,1)
    F(agent)=fobj(X(agent,:));
end
% Finds the agent with the best property code.
[Fbest,minp]=min(F);
Xbest=X(minp,:);
i=0;
% Main loop
for i=1:NumOfIteration 
    for agent=1:size(X,1)  
    % If the agent is the best agent, it will be rewarded with the bifurcated metamutation
        if agent==minp
            % Calculates absolute scattering Aavg and Xavg
            Xavg=mean(X);
            Ascat=mean(abs((X(agent,:)-Xavg)));
            % Calculates the reverse-weighted center
            Fc=abs(F);
            Xcenter=(1-Fc./(sum(Fc)))*X;
            % Bifurcated metamutation
            select=rand;
            if select <= 0.5
                Xn(agent,:)=Xcenter+(rand(1,size(X,2))-0.5)*Pout.* Ascat;
            else
                Xn(agent,:)=Xn(agent,:)+(rand(1,size(X,2))-0.5)*Pout.*Ascat;
            end
        else
            % Field-adapted differential crossover
            p1=rand(1,size(X,2)).*(Xbest-X(agent,:))*((abs(F(agent)-Fbest)))/(abs(F(agent))+abs(Fbest));
            % Field-aware mutation
            p2=(rand(1,size(X,2))-0.5)*(F(agent)/(F(agent)+Fbest));
            % Updates with a linear combination of both genetic processes.
            Xn(agent,:)=X(agent,:)+c1*p1+c2*p2;
        end
        % Check agents if they are within upper boundry (Ub) and lower boundray (Lb)  
        for var=1:size(X,2)
           if Xn(agent,var)>U(agent,var) || Xn(agent,var)<L(agent,var)
               % The properties, which are out of property ranges, are randomly placed into the interval(Lb,Ub).
               Xn(agent,:)=L(agent,:)+(U(agent,:)-L(agent,:))*rand;
           end
        end
    end
    
% Calculates the objective function values of new agents
        for agent=1:size(X,1)
            Fn(agent)=fobj(Xn(agent,:));
        end
% Form Old and new property old and new property code collection.
    Xg=[X;Xn];
    Fg=[F Fn];
% Natural selection is carried out according to minimum of objective function values.
    for agent=1:size(X,1)
           [Fmin ming]=min(Fg);
            X(agent,:)=Xg(ming,:);
            F(agent)=Fmin;
            Fg(ming)=Inf;
    end
    % Assigns the best agent of the current iteration
Fbest=F(1);
Xbest=X(1,:);    
BestAgentHist(i,:)=X(1,:);
BestCostHist(i)=F(1);

end
BestSol=BestAgentHist(i,:);
BestCost=BestCostHist(i);