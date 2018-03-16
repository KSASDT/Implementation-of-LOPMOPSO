clc;
clear;
close all;

%% Problem Definition

CostFunction=@(x) ZDT(x);      % Cost Function
m =2;
nVar=2;             % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin=0;          % Lower Bound of Variables
VarMax=1;          % Upper Bound of Variables

MaxVelocity = 0.2*(VarMax-VarMin);  % Upper bound on velocity
MinVelocity = -MaxVelocity;         % lower bound on velocity
	
	
%% MOPSO Parameters

MaxIt=200;           % Maximum Number of Iterations

nPop=200;            % Population Size

nRep=100;            % Archive Size

w=0.5;              % Inertia Weight
wdamp=0.99;         % Intertia Weight Damping Rate
c1=1;               % Personal Learning Coefficient
c2=2;               % Global Learning Coefficient

nGrid=7;            % Number of Grids per Dimension


beta=2;             % Leader Selection Pressure
gamma=2;            % Deletion Selection Pressure

mu=0.1;             % Mutation Rate

%% Initialization

empty_particle.Position=[];
empty_particle.Velocity=[];
empty_particle.Cost=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
empty_particle.IsDominated=[];
empty_particle.GridIndex=[];
empty_particle.GridSubIndex=[];

 % Create Population Array

pop=repmat(empty_particle,nPop,1);

 % Initialize Population Members
for i=1:nPop
    
	% Generate Random Solution
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
	 % Initialize Velocity
    pop(i).Velocity=zeros(VarSize);
    
	% Evaluation
    pop(i).Cost=CostFunction(pop(i).Position);
    
    
    % Update Personal Best
    pop(i).Best.Position=pop(i).Position;
    pop(i).Best.Cost=pop(i).Cost;
	
     
       % pdate Global Best
       % if particle(i).Best.Cost < GlobalBest.Cost
       %     GlobalBest = particle(i).Best;
       % end
    
end

% Determine Domination
pop=DetermineDomination(pop);
 
 
 % store nondominated solutions to rep
rep=pop(~[pop.IsDominated]);

Grid=CreateGrid(rep,nGrid);

for i=1:numel(rep)
    rep(i)=FindGridIndex(rep(i),Grid);
end


%set iteration counter=01

t=0

%% MOPSO Main Loop

while t < MaxIt

     % increment counter
     t=t+1;
	 
	 % set leader partical set rep2= null
	 
	 rep2 = []
	 
	  %calculating crowding distance
	 n = numel(rep);
	
	for i=1:numel(rep)
       pop(i).dist = 0;
    end
	 [n-2,2]= size(Y);
	 for k=1:m
       pop= sort(pop,k);
	   D(1)= D(n) = infinity
	   min=2;
	   for i=2: n-1
	     
         D(i) = pop(i+1)- pop(i-1);
		 If D(i)<D(i-1)
		    min=i;
       end
	 
	 Y = sort(D);
	 del(Y[n-2,1]);
	 D(min-1) = pop(min+1) - pop(min-2);
	 D(min+1) = pop(min+2) - (min-1);
	 
	 end

	 
function pop=DetermineDomination(pop)

    nPop=numel(pop);
    
    for i=1:nPop
        pop(i).IsDominated=false;
    end
    
    for i=1:nPop-1
        for j=i+1:nPop
            
            if Dominates(pop(i),pop(j))
               pop(j).IsDominated=true;
            end
            
            if Dominates(pop(j),pop(i))
               pop(i).IsDominated=true;
            end
            
        end
    end

end 
	 

function b=Dominates(x,y)

    if isstruct(x)
        x=x.Cost;
    end
    
    if isstruct(y)
        y=y.Cost;
    end

    b=all(x<=y) && any(x<y);

end

function z=ZDT(x)

    n=numel(x);

    f1=x(1);
    
    g=1+9/(n-1)*sum(x(2:end));
    
    h=1-sqrt(f1/g);
    
    f2=g*h;
    
    z=[f1
       f2];

end	 
	 
	 
 for i=1:nPop
        
		
		
        leader=SelectLeader(rep,beta);
		
		 % Update Velocity
        
        pop(i).Velocity = w*pop(i).Velocity ...
            +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position) ...
            +c2*rand(VarSize).*(leader.Position-pop(i).Position);
        
		
		% Apply Velocity Limits
            particle(i).Velocity = max(pop(i).Velocity, MinVelocity);
            particle(i).Velocity = min(pop(i).Velocity, MaxVelocity);
			
		 % Update Position
		
        pop(i).Position = pop(i).Position + pop(i).Velocity;
        
		% Apply Lower and Upper Bound Limits
        pop(i).Position = max(pop(i).Position, VarMin);
        pop(i).Position = min(pop(i).Position, VarMax);
        
		
		 % Evaluation
        pop(i).Cost = CostFunction(pop(i).Position);
        
       % Apply Mutation
      %  pm=(1-(it-1)/(MaxIt-1))^(1/mu);
       % if rand<pm
         %   NewSol.Position=Mutate(pop(i).Position,pm,VarMin,VarMax);
          %  NewSol.Cost=CostFunction(NewSol.Position);
           % if Dominates(NewSol,pop(i))
           %     pop(i).Position=NewSol.Position;
            %    pop(i).Cost=NewSol.Cost;

         %   elseif Dominates(pop(i),NewSol)
                % Do Nothing

          %  else
          %      if rand<0.5
           %         pop(i).Position=NewSol.Position;
           %%         pop(i).Cost=NewSol.Cost;
           %     end
       %     end
      %  end   
        
		
	 % Update Personal Best
	 
        if Dominates(pop(i),pop(i).Best)
            pop(i).Best.Position=pop(i).Position;
            pop(i).Best.Cost=pop(i).Cost;
            
        elseif Dominates(pop(i).Best,pop(i))
            % Do Nothing
            
        else
            if rand<0.5
                pop(i).Best.Position=pop(i).Position;
                pop(i).Best.Cost=pop(i).Cost;
            end
        end
        
    end
    
    % Add Non-Dominated Particles to REPOSITORY
    rep=[rep
         pop(~[pop.IsDominated])]; %#ok
    
    % Determine Domination of New Resository Members
    rep=DetermineDomination(rep);
    
    % Keep only Non-Dminated Memebrs in the Repository
    rep=rep(~[rep.IsDominated]);
    
    % Update Grid
    Grid=CreateGrid(rep,nGrid,alpha);

    % Update Grid Indices
    for i=1:numel(rep)
        rep(i)=FindGridIndex(rep(i),Grid);
    end
    
    % Check if Repository is Full
    if numel(rep)>nRep
        
        Extra=numel(rep)-nRep;
        for e=1:Extra
            rep=DeleteOneRepMemebr(rep,gamma);
        end
        
    end
    
    % Plot Costs
    figure(1);
    PlotCosts(pop,rep);
    pause(0.01);
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(numel(rep))]);
    
    % Damping Inertia Weight
    w=w*wdamp;
    
end

%% Resluts

