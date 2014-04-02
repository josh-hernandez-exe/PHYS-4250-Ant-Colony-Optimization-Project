function AntColony

par=setup();

World=SetupWorld(par);

%Global Best Solutions
GBTour = zeros(World.Enviro.n,1);
GBLength = inf ;

while (par.notDone)
    World.AntSystem=InitializeStart(World.AntSystem,World.Enviro.n,par);
    trip=2 ;
    
    while notDoneSearch(World.AntSystem,trip,par)   
        World = UpdateMemory(trip,World,par);
        trip=trip+1 ;       
    end % end while notDoneSearching
    
    World.AntSystem=MakeRoundTrip(World.AntSystem,par) ;
    
    World.AntSystem = CaculateTripLength(World,par);
    World.AntSystem.Rank=MakeRanking(World.AntSystem.Lengths,par);
    if par.GA.UseGA
        World.AntSystem=FindLocalMin(World,par);
    end
    World.Enviro = GloblePheromoneUpdate(World,par);
    
    par.T = 1 + par.T ;

	World.History.IBest(par.T)=World.AntSystem.Lengths(World.AntSystem.Rank(1));
    World.History.GBest(par.T)=GBLength;
    
    
    [GBTour,GBLength,par.show_prog]=GetResults(GBTour,GBLength,World,par);
    ShowBestLength(GBTour,GBLength,World.Enviro,par);
    par.notDone = AreWeDone(par);
    
end % end while notDone

% Show Final Results


figure('Name','Ant Ststem | Results','Numbertitle','off');
%Plot best route
subplot(2,2,1);
graph=World.Enviro.Nodes(GBTour,:);
if (World.Enviro.dim == 3)
    plot3(World.Enviro.Nodes(:,1),World.Enviro.Nodes(:,2),World.Enviro.Nodes(:,3),'b.',graph(:,1),graph(:,2),graph(:,3),'r.-');
else
    plot(World.Enviro.Nodes(:,1),World.Enviro.Nodes(:,2),'b.',graph(:,1),graph(:,2),'r.-');
end
title(sprintf('Total Distance = %1.4f, Iteration = %d',GBLength,par.T));
% Plot other results
subplot(2,2,2);
title(sprintf('Total Distance = %1.4f',GBLength));
plot(World.History.GBest,'b','LineWidth',2);
title('Global Best Solution History');
subplot(2,2,3);
plot(World.History.IBest,'r','LineWidth',2);
title('Iteration Best Solution History');
subplot(2,2,4);
imagesc(World.Enviro.Dist);
title('Distance Matrix');



end % end function AntColony
%******************************************************************
%******************************************************************
function World = SetupWorld(par)
Problem = SetupProblem();
World.Enviro = InitializeEnviro(Problem,par) ;
World.AntSystem = InitializeAntSystem(World.Enviro.n,par);
World.History.GBest=zeros(par.maxIteration,1);
World.History.IBest=zeros(par.maxIteration,1);
end % end function SetupWorld
%******************************************************************
%******************************************************************
function param = InitializeAntSystem(n,par)
m=par.AntNum;

AntTours = zeros(m,n+1);
ToursLength = zeros(m,1);
Move = ones(m,1);
Rank = zeros(m,1);
param = struct('Tours',AntTours,'Lengths',ToursLength,'Move',Move,'Rank',Rank);
end % function InitializeAntSystem
%******************************************************************
%******************************************************************
function param = InitializeEnviro(Problem,par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fix after as it depends on problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.n=Problem.n;
param.dim=Problem.dim;
param.Nodes=Problem.map;
param.maxDist=max(max(Problem.dist));
%param.Dist=Problem.dist ./ param.maxDist; % Normalized distance
param.Dist=Problem.dist;
param.Tau=(ones(param.n,param.n)-eye(param.n,param.n))*par.tau_0;

% for i=1:param.n
%     for j=1:param.n
%         param.Tau(i,j)=(i-1)*param.n+j;
%     end
% end

param.ShouldAntGo = @(AntSystem,n,trip,par) Problem.ShouldAntGo(AntSystem,n,trip,par);
end % end InitializeEnviro
%******************************************************************
%******************************************************************
function isNotDone = notDoneSearch(AntSystem,trip,par)
isNotDone=1 ;

if ( all(AntSystem.Move(:,:)==0) || trip > par.maxTrips )
    isNotDone=0 ;
end % end if
    
end % end function constrantFunction
%******************************************************************
%******************************************************************
function param = InitializeStart(AntSystem,n,par)
m=par.AntNum;

AntSystem.Tours=zeros(m,n+1);
AntSystem.Tours(:,1) = ceil(n*rand(m,1));
AntSystem.Lengths = zeros(m,1);
AntSystem.Move = ones(m,1);
AntSystem.Rank = zeros(m,1);
param=AntSystem ;
end % end function InitStartPoint
%******************************************************************
%******************************************************************
function P = UpdateProbability(trip,World,par)

CurrentPoint = World.AntSystem.Tours(:,trip-1);
VisitedPoints=zeros(par.AntNum,trip-1);
tau=zeros(par.AntNum,World.Enviro.n);
dist=zeros(par.AntNum,World.Enviro.n);

%This calculates the distance to all possible nodes
%using the past knowledge of the tour

for ant=1:par.AntNum
    CP=CurrentPoint(ant);
    VisitedPoints(ant,:) = World.AntSystem.Tours(ant, 1:trip-1);
    tau(ant,:) = World.Enviro.Tau(CP,:);
    tau(ant,VisitedPoints(ant,:))=0;
    dist(ant,:) = World.Enviro.Dist(CP,:);
    dist(ant,CP)=1;
end

P = ProbabilityFunction(tau,dist,par);

P_sum = sum(P'); %note that sum command sums vertically

NotZero= (P_sum~= 0);

for ant=1:par.AntNum
    if NotZero(ant)
        P(ant,:) = P(ant,:)./P_sum(ant);
    else 
        NoVisitedNodes = setdiff(1:World.Enviro.n,VisitedPoints(ant,:));
        P(ant,NoVisitedNodes) = 1/length(NoVisitedNodes); 
    end% end if else
end
end % end function UpdateProbability
%******************************************************************
%******************************************************************
function P = ProbabilityFunction(tau,dist,par)
P = (tau.^par.alpha).*((1./dist).^par.beta);
end % function ProbabilityFunction
%******************************************************************
%******************************************************************
function Select = ChoiceFunction(P,par)
%Calculates which node to visit
[num,m] = size(P);
Select = zeros(1,num);
P_sum = sum(P');
r = rand(1,num);
for ant=1:num
    if rand<=par.q_0
         WhereIsMax=find(max(P(ant,:))==P(ant,:));
         if length(WhereIsMax)>1
             Select(ant)=WhereIsMax(ceil(length(WhereIsMax)*rand));
         elseif length(WhereIsMax)==1
             Select(ant)=WhereIsMax;
         else
             Select(ant)=ceil(num*rand);
         end
    else
        flag = (1-P_sum(ant)<=1e-5);
        sumP = 0;
        j = ceil(m*rand); 
        while (sumP<r(ant)) && flag
            sumP = sumP + P(ant,mod(j-1,m)+1);
            j = j+1;
        end % end while
        Select(ant) = mod(j-2,m)+1;
    end % if else
end % end for
end % end function ChoiceFunction
%******************************************************************
%******************************************************************
function param = UpdateMemory(trip,World,par)

%         for ant = 1 : par.AntNum          
%             if(World.AntSystem.Move(ant))     
%                 World.AntSystem = UpdateMemory(trip,ant,World,par);
%             end % end if Move
%         end % end for ant

        
P = UpdateProbability(trip,World,par);
NextPoint = ChoiceFunction(P,par);
World.AntSystem.Tours(:,trip) = NextPoint;
World.Enviro=Online_step_by_step_pheromone(World,trip,par);
World.AntSystem.Move=World.Enviro.ShouldAntGo(World.AntSystem,World.Enviro.n,trip,par);
param = World ;
end
%******************************************************************
%******************************************************************
function Enviro=Online_step_by_step_pheromone(World,trip,par)
% alters the pheromone during the creation of solutinos
if par.Reduction
    
    m=par.AntNum;
    
    for ant=1:m
        pointA=World.AntSystem.Tours(ant,(trip-1));
        pointB=World.AntSystem.Tours(ant,trip);
        World.Enviro.Tau(pointA,pointB) = (1-par.phi)*World.Enviro.Tau(pointA,pointB) + par.phi*par.tau_0 ;
    end
end % end if

Enviro=World.Enviro;
end
%******************************************************************
%******************************************************************
function AntSystem = CaculateTripLength(World,par)

Lengths = zeros(par.AntNum,1);
for k=1:par.AntNum
    for i=1:(length(World.AntSystem.Tours(k,:))-1)
        pointA=World.AntSystem.Tours(k,i);
        pointB=World.AntSystem.Tours(k,i+1);
        Lengths(k)=Lengths(k)+ World.Enviro.Dist(pointA,pointB);
    end
end
World.AntSystem.Lengths = Lengths;

AntSystem = World.AntSystem ;

end % end function CaculateTripLength
%******************************************************************
%******************************************************************
function Rank = MakeRanking(Lengths,par)

Array = [Lengths,(1:length(Lengths))'];

Array = MergeSort(Array);

Rank = Array(:,2);
end %
%******************************************************************
%******************************************************************
%**************************MERGE SORT CODE*************************
%******************************************************************
function result = MergeSort(Array)
size=length(Array(:,1));
if size<=1
    result = Array;
else
    middle = floor(size/2);
    left = Array(1:middle,:);
    right = Array(middle+1:size,:);
    left = MergeSort(left);
    right = MergeSort(right);
    result = merge(left,right);
end % else if
end %
%******************************************************************
%******************************************************************
function result = merge(left,right)
sizeL=length(left(:,1));
sizeR=length(right(:,1));
result = zeros(sizeL+sizeR,2);
counter=1;
while (sizeL > 0 || sizeR > 0 )
    if (sizeL > 0 && sizeR > 0 )
        if  (left(1,1) <= right(1,1) )
            
            result(counter,:)= left(1,:);
            left =left(2:sizeL,:); %The rest of left
        else
            result(counter,:)= right(1,:);        
            right =right(2:sizeR,:); %The rest of right
        end % end if
    elseif sizeL > 0
        result(counter,:)= left(1,:);
        left =left(2:sizeL,:); %The rest of left
    elseif sizeR > 0
        result(counter,:)= right(1,:);        
        right =right(2:sizeR,:); %The rest of right
    end %end if
    sizeL=length(left(:,1));
    sizeR=length(right(:,1));
    counter=1+counter;
end % end while
end % end function merge
%******************************************************************
%******************************************************************
%***********************END MERGE SORT CODE************************
%******************************************************************
function AntSystem=FindLocalMin(World,par)
%Use Genetic Algorithm to find local min of Best Solution
BestAnts=World.AntSystem.Rank;


num_iter =par.GA.num_iter ;
n=World.Enviro.n;
Dist=World.Enviro.Dist;

% Sanity Checks
pop_size = 4*par.GA.mod4*par.GA.pickBest;


% Initialize the Population
pop = zeros(pop_size,n+1);
tmp_pop = zeros(4,n+1);
new_pop = zeros(pop_size,n+1);
total_dist=zeros(1,pop_size);


for i = 1:par.GA.pickBest
    for j=1:4*par.GA.mod4
        pop(4*par.GA.mod4*(i-1)+j,:) = World.AntSystem.Tours(BestAnts(i),:) ;
        total_dist(4*par.GA.mod4*(i-1)+j)=World.AntSystem.Lengths(BestAnts(i),:);
    end
end

% Run the GA
for iter = 1:num_iter   
    %Algorithm Operators
    rand_pair = randperm(pop_size);
    for p = 4:4:pop_size
        rtes = pop(rand_pair(p-3:p),:);
        dists = total_dist(rand_pair(p-3:p));
        [ignore,idx] = min(dists);
        best_of_4_rte = rtes(idx,:);
        ins_pts = sort(ceil(n*rand(1,2)));
        I = ins_pts(1);
        J = ins_pts(2);
        for k = 1:4 % Mutate the Best to get Three New Routes
            tmp_pop(k,:) = best_of_4_rte;
            switch k
                case 2 % Flip
                    tmp_pop(k,I:J) = fliplr(tmp_pop(k,I:J));
                case 3 % Swap
                    tmp_pop(k,[I J]) = tmp_pop(k,[J I]);
                case 4 % Slide
                    tmp_pop(k,I:J) = tmp_pop(k,[I+1:J I]);
                otherwise % Do Nothing
            end
        end
        new_pop(p-3:p,:) = tmp_pop;
    end
    pop = new_pop; 
    % Evaluate Each Population Member (Calculate Total Distance)
    for p = 1:pop_size
        d = Dist(pop(p,n+1),pop(p,1)); % Closed Path
        for k = 1:n
            d = d + Dist(pop(p,k),pop(p,k+1));
        end
        total_dist(p) = d;
    end    
end

Rank=MakeRanking(total_dist',par);

for ant=1:par.GA.pickBest
    World.AntSystem.Tours(BestAnts(par.AntNum+1-ant),:)=pop(Rank(ant),:);
    World.AntSystem.Lengths(BestAnts(par.AntNum+1-ant))=total_dist(Rank(ant));
end %end for


World.AntSystem.Rank=MakeRanking(World.AntSystem.Lengths,par);

AntSystem = World.AntSystem ;

end
%******************************************************************
%******************************************************************
function Enviro = GloblePheromoneUpdate(World,par)
m=par.AntNum;
n=World.Enviro.n;
rho=par.rho;

Tours = World.AntSystem.Tours;
Lengths = World.AntSystem.Lengths;
deltaTau= (1)./Lengths ;
sumDeltaTau_reg = zeros(n,n);
sumDeltaTau_rank = zeros(n,n);

%These are based of the forumlas that determine the alteration
%of pheromone levels
if par.RegPheromone
    for k=1:m
        for i=1:n
            sumDeltaTau_reg(Tours(k,i),Tours(k,i+1))=sumDeltaTau_reg(Tours(k,i),Tours(k,i+1))+(par.Q)*deltaTau(k);
        end % end for
    end % end for
end % end if

if par.RankPheromone
    if (par.sigma>m)
        par.sigma=m;
    end % end if
    for k=1:(par.sigma-1)
        ant_k=World.AntSystem.Rank(k);
        for i=1:n
            pointA=Tours(ant_k,i);
            pointB=Tours(ant_k,i+1);
            sumDeltaTau_rank(pointA,pointB)=sumDeltaTau_rank(pointA,pointB)+(par.sigma-k)*(par.QRank)*deltaTau(ant_k);
            if(k==1)
                sumDeltaTau_rank(pointA,pointB)=sumDeltaTau_rank(pointA,pointB)+par.QBest*par.sigma*deltaTau(ant_k);
            end % end if
        end % end for
    end % end for
end % end if

World.Enviro.Tau = (1-rho)*World.Enviro.Tau + sumDeltaTau_reg + sumDeltaTau_rank ;

if (par.tau_max <inf)
    OverMax=World.Enviro.Tau>par.tau_max;
    if any(OverMax)
        World.Enviro.Tau(OverMax)=par.tau_max;
    end
end

if (par.tau_min >0)
    UnderMin= World.Enviro.Tau<par.tau_min;
    UnderMin = logical(UnderMin) ;
    World.Enviro.Tau(UnderMin)=par.tau_min;
end
Enviro = World.Enviro ;

% if (mod(par.T,par.DispInterval)==0)
%     disp('Biggest Tau');
%     LessThan = World.Enviro.Tau<=par.tau_max;
%     disp(max(max(World.Enviro.Tau(LessThan))));
%     disp('Biggest Tau < MaxTau');
%     LessThan = World.Enviro.Tau<par.tau_max;
%     disp(max(max(World.Enviro.Tau(LessThan))));
%     disp('Smallest Tau');
%     MoreThan=World.Enviro.Tau>=par.tau_min;
%     disp(min(min(World.Enviro.Tau(MoreThan))));
%     disp('Smallest Tau > MinTau');
%     MoreThan=World.Enviro.Tau>par.tau_min;
%     disp(min(min(World.Enviro.Tau(MoreThan))));
% end


end % end function GloblePheromoneUpdate
%******************************************************************
%******************************************************************
function notDone = AreWeDone(par)
notDone = 1 ;
if (par.T>=par.maxIteration )
    notDone = 0;
end
end % end function AreWeDone
%******************************************************************
%******************************************************************
function param = MakeRoundTrip(AntSystem,par)
if par.MakeRoundTrip
    
    %LastTrip=zeros(par.AntNum,1)
    %for i=1:par.AntNum
    %    LastTrip(i)=length(AntSystem.Tours(i,:));
    %end % end for
    
    n=length(AntSystem.Tours(1,:));
    
    AntSystem.Tours(:,n) = AntSystem.Tours(:,1);
end % end if
param=AntSystem ;
end
%******************************************************************
%******************************************************************
function [GBTour,GBLength,show_prog]=GetResults(PBTour,PBLength,World,par)

BestAnt= World.AntSystem.Rank(1);
IBTour = World.AntSystem.Tours(BestAnt,:);
IBLength = World.AntSystem.Lengths(BestAnt,:);

if (IBLength<=PBLength)
    show_prog=1;
    GBTour =IBTour;
    GBLength = IBLength ;
else
    show_prog=0;
    GBTour = PBTour;
    GBLength = PBLength;
end % end if elseif else

end %end function GetResults
%******************************************************************
%******************************************************************
function ShowBestLength(GBTour,GBLength,Enviro,par)

if par.show_prog %&& (mod(par.T,par.DispInterval)==0)
    %Plot the Best Route
    %GBLength=GBLength*Enviro.maxDist;
    figure(par.fig);
    graph=Enviro.Nodes(GBTour,:);
    if (Enviro.dim == 3)
        plot3(Enviro.Nodes(:,1),Enviro.Nodes(:,2),Enviro.Nodes(:,3),'b.',graph(:,1),graph(:,2),graph(:,3),'r.-');
    else
        plot(Enviro.Nodes(:,1),Enviro.Nodes(:,2),'b.',graph(:,1),graph(:,2),'r.-');
    end
    title(sprintf('Total Distance = %1.4f, Iteration = %d',GBLength,par.T));

    pause(0.01);

end


end
%****************************************************************
%*******************************SETUP****************************
%****************************************************************
function par=setup()

par = OptionSetup();

par.T = 0 ;
par.notDone= 1 ; % termination criterion
end % function setup
%******************************************************************
%******************************************************************
function par=OptionSetup()
%******************************************************************
%This is an outline on how to change the setup to match
% a particular implamentation of the Ant System
%
% For Ant System
%
%         par.tau_min=0;
%         par.tau_max=inf;
%         par.tau_0=1e-3;
% 
%         par.RegPheromone =1 ;
%         par.RankPheromone = 0;
%         par.Reduction=0;
%         par.GA.UseGA=0;
% 
% For Max Min Ant System
%
%         par.tau_min=0.1;  % they can be anything really
%         par.tau_max=5;    % as long as min < max
%         par.tau_0=par.tau_max;
% 
%         par.RegPheromone =0 ;
%         par.RankPheromone = 1;
%         par.Reduction=0;
%         par.GA.UseGA=0;
% 
%         par.sigma = 1;
% 
% For Ant System Rank
%
%         par.tau_min=0;
%         par.tau_max=inf;
%         par.tau_0=1e-3;
% 
%         par.RegPheromone =0 ;
%         par.RankPheromone = 1;
%         par.Reduction=0;
%         par.GA.UseGA=0;
% 
%         par.sigma = 5;  % can be anything as long as par.sigma < par.AntNum
% 
% For Ant Colony System
%
%         par.tau_min=0;
%         par.tau_max=inf;
%         par.tau_0=1e-3;
% 
%         par.RegPheromone =0 ;
%         par.RankPheromone = 1;
%         par.Reduction=1;
%         par.GA.UseGA=1;
% 
%         par.sigma = 1;
% 
%         par.GA.pickBest = 1;
%         par.GA.mod4=1;
%         par.GA.num_iter = 1;
% 
%For Hybrid System
%
%         par.tau_min=0;
%         par.tau_max=inf;
%         par.tau_0=1e-3;
% 
%         par.RegPheromone =0 ;
%         par.RankPheromone = 1;
%         par.Reduction=1;
%         par.GA.UseGA=1;
% 
%         par.sigma = 5;
% 
%         par.GA.pickBest = 3;
%         par.GA.mod4=3;
%         par.GA.num_iter = 3;
 


par.maxIteration = 2000 ;
par.maxTrips=200; % Just in case the creation of a solution excedes the regular amount
par.AntNum= 48 ; % Should be equal to the dimension
par.DispInterval = 10;

par.show_prog=1;
par.fig=figure('Name','Ant System | Current Best Solution','Numbertitle','off');

% Affectes Pheromone
par.alpha = 1;
par.beta = 5;
par.rho = 0.5; % rho \in (0,1]
par.Q=10;
par.QRank=10;
par.QBest=1;
par.q_0=0.15;  % q_0 \in (0,1] 

par.tau_min=0;    % the min tau
par.tau_max=inf;  % the max tau
par.tau_0=1e-3;   % the initial tau

% ON / OFF Parameters
par.MakeRoundTrip=1 ;
par.RegPheromone =0 ;    % use regular pheromone update rules
par.RankPheromone = 1;   % use rank pheromone update rules
par.Reduction=1;         % use online pheromone update rules
par.GA.UseGA=1;          % use GA as local optimizer

% used when  RankPheromone==1
par.sigma = 5;      % the top number of ants that will be picked for ranking

% used when Reduction=1
par.phi=0.25; % phi \in (0,1]

%Genetic Algorithm Parameters
par.GA.pickBest = 3;  % the top number of ants that will be locally optimized
par.GA.mod4=3;        % the number of copies of each solution will be added to poulation
par.GA.num_iter = 3;  % how many iterations the GA will run for

end % function OptionSetup
%******************************************************************
%******************************************************************
%******************PROBLEM SPECIFIC FUNTIIONS *********************
%******************************************************************
%******************************************************************
function prob=SetupProblem()
prob.map = MapSetup();
prob.dist=distanceFuntion(prob.map);


[prob.n,prob.dim]=size(prob.map);

prob.MaxCost = 0 ;

prob.ShouldAntGo = @(AntSystem,n,trip,par) ShouldAntGo(AntSystem,n,trip,par);

end % end probSetup
%******************************************************************
%******************************************************************
function map = MapSetup()
% This really defines the problem
%Odyssey of Ulysses with 48 points
map(48,1:2)=[0,0];

map(1,1:2) = [ 6734, 1453];
map(2,1:2) = [ 2233, 10	 ];
map(3,1:2) = [ 5530, 1424];
map(4,1:2) = [ 401 , 841 ];
map(5,1:2) = [ 3082, 1644];
map(6,1:2) = [ 7608, 4458];
map(7,1:2) = [ 7573, 3716];
map(8,1:2) = [ 7265, 1268];
map(9,1:2) = [ 6898, 1885];
map(10,1:2)= [ 1112, 2049];
map(11,1:2)= [ 5468, 2606];
map(12,1:2)= [ 5989, 2873];
map(13,1:2)= [ 4706, 2674];
map(14,1:2)= [ 4612, 2035];
map(15,1:2)= [ 6347, 2683];
map(16,1:2)= [ 6107, 669 ];
map(17,1:2)= [ 7611, 5184];
map(18,1:2)= [ 7462, 3590];
map(19,1:2)= [ 7732, 4723];
map(20,1:2)= [ 5900, 3561];
map(21,1:2)= [ 4483, 3369];
map(22,1:2)= [ 6101, 1110];
map(23,1:2)= [ 5199, 2182];
map(24,1:2)= [ 1633, 2809];
map(25,1:2)= [ 4307, 2322];
map(26,1:2)= [ 675 , 1006];
map(27,1:2)= [ 7555, 4819];
map(28,1:2)= [ 7541, 3981];
map(29,1:2)= [ 3177, 756 ];
map(30,1:2)= [ 7352, 4506];
map(31,1:2)= [ 7545, 2801];
map(32,1:2)= [ 3245, 3305];
map(33,1:2)= [ 6426, 3173];
map(34,1:2)= [ 4608, 1198];
map(35,1:2)= [ 23.0, 2216];
map(36,1:2)= [ 7248, 3779];
map(37,1:2)= [ 7762, 4595];
map(38,1:2)= [ 7392, 2244];
map(39,1:2)= [ 3484, 2829];
map(40,1:2)= [ 6271, 2135];
map(41,1:2)= [ 4985, 140 ];
map(42,1:2)= [ 1916, 1569];
map(43,1:2)= [ 7280, 4899];
map(44,1:2)= [ 7509, 3239];
map(45,1:2)= [ 10.0, 2676];
map(46,1:2)= [ 6807, 2993];
map(47,1:2)= [ 5185, 3258];
map(48,1:2)= [ 3023, 1942];
end % function MapSetup
%******************************************************************
%******************************************************************
function D=distanceFuntion(map)
% One can put in there own distance matrix if they like
x = map(:,1);
y= map(:,2);
n=length(map(:,1)); % map(:,2) would be the same

Dx =  ones(n,1) * x' - x * ones(1,n);
Dy =  ones(n,1) * y' - y * ones(1,n);

Dsq = Dx .* Dx + Dy .* Dy;
D = sqrt(Dsq);
end % end function dist
%******************************************************************
%******************************************************************
function shouldGo = ShouldAntGo(AntSystem,n,trip,par)
% this is technically problem specific 
% and will drasticly change depending on the problem
shouldGo= (AntSystem.Tours(:,n,:)== 0);

end % end Should AntGo