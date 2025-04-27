function [Xsun,Gsun,convergence]=POA(N,T_max,low,up,dim,fobj)
% Define main factors
G=1;      % The gravitational constant
a=2;      % The constant to calculate the mass
b=5;      % The constant parameter is generally chosen from 2 to 10
R0=1000;  % The factor for dividing the search space

% Initialize Gsun, Xsun, MF
Gsun=inf;
Xsun=Gsun*ones(1,dim);
fitness=Gsun*ones(N,1);
MF=zeros(1,N); % Moment Force

% Initialize the positions
X=low+(up-low).*rand(N,dim);
X_new=zeros(N,dim+1);


% Define Rmin value
Rmin=sqrt(sum((up-low).^2))/R0; % As shown in Eq.(6)

% Determine inital Gsun and Xsun
for i=1:N
    fitness(i)=fobj(X(i,:));
end
X_new= sortrows([X fitness],dim+1);
Xsun=X_new(1,1:dim);
Gsun=X_new(1,dim+1);
fitmax=(max(fitness));
fitmin=(min(fitness));
alpha=abs(fitmax-fitmin);

% Main loop
for t=1:T_max
    Ms=1/a^(Gsun/alpha); % The mass of the Sun
    for i=1:N
        R(i)=POAdistance(X(i,:),Xsun);  % Distance between the Sun and Planet(i) as shown in Eq.(3)
        if R(i)<Rmin
            MF(i)=0;
        else
            Mass=1/a^(fitness(i)/alpha); % The mass of planet(i) 
            MF(i)=G*Mass*Ms/R(i); % Calculate Moment Force as shown in Eq.(2)
        end
    end
    
    beta=MF/max(MF); % beta = MF(i)/MF(max) 
    
    % Update c: a decrease factor from 2 to 1
    c=2-t/T_max;
    
    for i=1:N
        r1=rand(1,dim); % a random number in (0;1)
        if R(i) > Rmin
            X(i, :) = X(i, :) +b*beta(i)*r1.*(Xsun-X(i, :)); % Update location of Planet(i) as shown in Eq.(4)
        else
            r2=0.2*randn+0.5; % r2 is a normal distribution function
            X(i, :) = X(i, :) +c*r1.*(r2*Xsun - X(i, :)); % Update location of Planet(i) as shown in Eq.(5)
        end
        
        % Boudary of the Problem
        xmin=X(i,:)<low;
        xmax=X(i,:)>up;
        X(i,:)=xmin.*low+xmax.*up+(ones(1,dim)-xmin-xmax).*X(i,:);

        %Evaluate new solutions
        fitness(i)=fobj(X(i,:));
    end
    
    % Update the Sun and parameters to calculate the mass
    X_new= sortrows([X fitness],dim+1);
    if X_new(1,dim+1)<Gsun
        Xsun=X_new(1,1:dim);
        Gsun=X_new(1,dim+1);
    end     
    fitmax=(max(fitness));
    fitmin=(min(fitness));
    alpha=abs(fitmax-fitmin);
    
    convergence(t)=(Gsun);

end
end

function dist = POAdistance(x1, x2)
    % Calculate Euclidean distance between two points
    dist = sqrt(sum((x1 - x2).^2));
end