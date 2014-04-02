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
        for ant = 1 : par.AntNum          
            if(World.AntSystem.Move(ant))     
                World.AntSystem = UpdateMemory(trip,ant,World,par);
            end % end if Move
        end % end for ant
        trip=trip+1 ;       
    end % end while notDoneSearching
    
    World.AntSystem=MakeRoundTrip(World.AntSystem,par) ;
    
    World.AntSystem = CaculateTripLength(World,par);
    World.AntSystem.Rank=MakeRanking(World.AntSystem,par);
    World.Enviro = GloblePheromoneUpdate(World,par);
    
    par.T = 1 + par.T ;
    
    [GBTour,GBLength]=GetResults(GBTour,GBLength,World,par);
    ShowBestLength(GBTour,GBLength,World.Enviro.Nodes,par);
    
    par.notDone = AreWeDone(par);
    
end % end while notDone
end % end function AntColony
%******************************************************************
%******************************************************************
function World = SetupWorld(par)
Problem = SetupProblem();
World.AntSystem = InitializeAntSystem(Problem.nPoints,par);
World.Enviro = InitializeEnviro(Problem,par) ;
end % end function SetupWorld
%******************************************************************
%******************************************************************
function param = InitializeAntSystem(n,par)
m=par.AntNum;

if par.MakeRoundTrip
    AntTours = zeros(m,n+1);
else
    AntTours = zeros(m,n);
end
ToursLength = zeros(m,1);
Move = ones(m,1);
Rank = zeros(m,1);
param = struct('Tours',AntTours,'Lengths',ToursLength,'Move',Move,'Rank',Rank);
end % function InitializeAntSystem
%******************************************************************
%******************************************************************
function param = InitializeEnviro(Problem,par)
n=Problem.nPoints;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fix after as it depends on problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MatrixTau = (ones(n,n)-eye(n,n))*par.tau_0;
Dist= Problem.dist;
Cost=Problem.cost;
param=struct('Nodes',Problem.map,'Dist',Dist,'Tau',MatrixTau,'Cost',Cost);
param.n=n;
param.ShouldAntGo = @(ant,AntSystem,n,par) Problem.ShouldAntGo(ant,AntSystem,n,par);
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

if par.MakeRoundTrip
    AntSystem.Tours=zeros(m,n+1);
else
    AntSystem.Tours=zeros(m,n);
end

AntSystem.Tours(:,1) = randint(m,1,[1,n]);
AntSystem.Lengths = zeros(m,1);
AntSystem.Move = ones(m,1);
AntSystem.Rank = zeros(m,1);
param=AntSystem ;
end % end function InitStartPoint
%******************************************************************
%******************************************************************
function P = UpdateProbability(trip,ant,World,par)

CurrentPoint = World.AntSystem.Tours(ant, trip-1);
VisitedPoints = World.AntSystem.Tours(ant, 1:trip-1);
tau = World.Enviro.Tau(CurrentPoint,:);
tau(1,VisitedPoints) = 0;
dist = World.Enviro.Dist(CurrentPoint,:)*World.Enviro.Cost;
dist(1,CurrentPoint) = 1;

P = ProbabilityFunction(tau,dist,par);

if sum(P) ~= 0
    P = P/sum(P);
else 
    NoVisitedNodes = setdiff(1:World.Enviro.n,VisitedPoints);
    P(1,NoVisitedNodes) = 1/length(NoVisitedNodes); 
end% end if else
end % end function UpdateProbability
%******************************************************************
%******************************************************************
function P = ProbabilityFunction(tau,dist,par)
P = (tau.^par.alpha).*((1./dist).^par.beta);
end % function ProbabilityFunction
%******************************************************************
%******************************************************************
function Select = ChoiceFunction(P,num)
m = length(P);
flag = (1-sum(P)<=1e-5);
Select = zeros(1,num);

r = rand(1,num);
for i=1:num
    sumP = 0;
    j = ceil(m*rand); 
    while (sumP<r(i)) && flag
        sumP = sumP + P(mod(j-1,m)+1);
        j = j+1;
    end % end while
    Select(i) = mod(j-2,m)+1;
end % end for

end % end function ChoiceFunction
%******************************************************************
%******************************************************************
function param = UpdateMemory(trip,ant,World,par)
P = UpdateProbability(trip,ant,World,par);
NextPoint = ChoiceFunction(P,1);
World.AntSystem.Tours(ant,trip) = NextPoint;
World.AntSystem.Move(ant)=World.Enviro.ShouldAntGo(ant,World.AntSystem,World.Enviro.n,par);
param = World.AntSystem ;
end
%******************************************************************
%******************************************************************
function AntSystem = CaculateTripLength(World,par)

Lengths = zeros(par.AntNum,1);
for k=1:par.AntNum
    for i=1:(length(World.AntSystem.Tours(k,:))-1)
        for j=1:(length(World.AntSystem.Tours(k,:))-1)
            if (i~=j)
                pointA=World.AntSystem.Tours(k,i);
                pointB=World.AntSystem.Tours(k,j);
                Lengths(k)=Lengths(k)+World.Enviro.Cost(i,j)*World.Enviro.Dist(pointA,pointB);
            end
        end
    end
end
World.AntSystem.Lengths = Lengths;

AntSystem = World.AntSystem ;

end % end function CaculateTripLength
%******************************************************************
%******************************************************************
function Rank = MakeRanking(AntSystem,par)

Array = [AntSystem.Lengths,(1:par.AntNum)'];

Array = MergeSort(Array);

Rank = Array(:,2);
end %
%******************************************************************
%******************************************************************
%**************************MERGE SORT CODE*************************
%******************************************************************
function result = MergeSort(Array)
size=length(Array(:,1));
middle = floor(size/2) ;

if size<=1
    result = Array;
else
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
    sizeL=length(left(:,1));
    sizeR=length(right(:,1));
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
    counter=1+counter;
end % end while
end % end function merge
%******************************************************************
%******************************************************************
%***********************END MERGE SORT CODE************************
%******************************************************************
function Enviro = GloblePheromoneUpdate(World,par)
m=par.AntNum;
n=World.Enviro.n;
rho=par.rho;

Tours = World.AntSystem.Tours;
Lengths = World.AntSystem.Lengths;

if par.RegPheromone
    deltaTau= (par.Q)./Lengths ;
    sumDeltaTau = zeros(n,n);
    for k=1:m
        for i=1:length(Tours(k,:))-1
            sumDeltaTau(Tours(k,i),Tours(k,i+1))=sumDeltaTau(Tours(k,i),Tours(k,i+1))+deltaTau(k);
        end % end for
    end % end for
    World.Enviro.Tau = (1-rho)*World.Enviro.Tau + sumDeltaTau ;
end % end if

if par.RankPheromone
    deltaTau= (par.QRank)./Lengths ;
	sumDeltaTau = zeros(n,n);
    
    if (par.sigma>m)
        par.sigma=m;
    end % end if
    for k=1:(par.sigma-1)
        ant_k=World.AntSystem.Rank(k);
        for i=1:length(Tours(ant_k,:))-1
            pointA=Tours(ant_k,i);
            pointB=Tours(ant_k,i+1);
            sumDeltaTau(pointA,pointB)=sumDeltaTau(pointA,pointB)+(par.sigma-k)*deltaTau(ant_k);
            if(k==1)
                sumDeltaTau(pointA,pointB)=sumDeltaTau(pointA,pointB)+par.sigma*deltaTau(ant_k);
            end % end if
        end % end for
    end % end for
    World.Enviro.Tau = (1-rho)*World.Enviro.Tau + sumDeltaTau ;
end % end if


if (par.tau_max >0)
    OverMax=World.Enviro.Tau>par.tau_max;
    World.Enviro.Tau(OverMax)=par.tau_max;
end

if (par.tau_min >0)
    UnderMin= World.Enviro.Tau<par.tau_min;
    World.Enviro.Tau(UnderMin)=par.tau_min;
end
Enviro = World.Enviro ;

% if (mod(par.T,par.DispInterval)==0)
% 	disp('Biggest Tau < MaxTau');
%     LessThan = World.Enviro.Tau<par.tau_max;
%     disp(max(max(World.Enviro.Tau(LessThan))));
% 	disp('Smallest Tau > MinTau');
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
function [GBTour,GBLength]=GetResults(PBTour,PBLength,World,par)

BestAnt= World.AntSystem.Rank(1);
IBTour = World.AntSystem.Tours(BestAnt,:);
IBLength = World.AntSystem.Lengths(BestAnt,:);

if (IBLength<=PBLength)
    GBTour =IBTour;
    GBLength = IBLength ;
else
    GBTour = PBTour;
    GBLength = PBLength;
end % end if elseif else

end %end function GetResults
%******************************************************************
%******************************************************************
function ShowBestLength(GBTour,GBLength,Nodes,par)
if (mod(par.T,par.DispInterval)==0)
    disp({'Iteration :',par.T});
    disp('Facility');
    disp(Nodes(1,:));
    disp('Location');
    disp(GBTour);
    disp({'Lowest Cost',GBLength});
end % end if

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

par.maxIteration = 10000 ;
par.maxTrips=100;
par.AntNum= 64 ; % Should be equal to the dimension
par.lambda = 0.15;
par.DispInterval = 10;

% Affectes Pheromone
par.alpha = 1 ;
par.beta = 5;
par.rho = 0.6 ;
par.Q=1;
par.QRank=2;

par.tau_min=5;
par.tau_max=20;
par.tau_0=20;

% ON / OFF Parameters
par.MakeRoundTrip=0 ;
par.RegPheromone =1 ;
par.RankPheromone =1;


% used when  RankPheromone==1
par.sigma = 100;

end % function OptionSetup
%******************************************************************
%******************************************************************
%******************PROBLEM SPECIFIC FUNTIIONS *********************
%******************************************************************
%******************************************************************
function prob=SetupProblem()
prob.map=MapSetup();
prob.dist=distanceFuntion();
prob.cost=costFunction();

prob.MaxCost = 0 ;

prob.nPoints=length(prob.dist(:,1));
% map(:,2) would give the same result
% both should have the same dimension

prob.ShouldAntGo = @(ant,AntSystem,n,par) ShouldAntGo(ant,AntSystem,n,par);

end % end probSetup
%******************************************************************
%******************************************************************
function map=MapSetup()

map(1:8)=[1:8];
end
%******************************************************************
%******************************************************************
function D=distanceFuntion()
D = [0 1 2 3 1 2 3 4;
      5 0 1 2 2 1 2 3;
      2 3 0 1 3 2 1 2;
      4 0 0 0 4 3 2 1;
      1 2 0 5 0 1 2 3;
      0 2 0 2 10 0 1 2;
      0 2 0 2 0 5 0 1;
      6 0 5 10 0 1 10 0];
size=length(D(1,:));
  
for i=1:size
    for j=1:size
        if (i<j)       
            D(i,j)=D(j,i);
        end % end if
    end %end for
end % end for

end % end function dist
%******************************************************************
%******************************************************************
function D=costFunction()
D = [0 1 2 3 1 2 3 4;
      5 0 1 2 2 1 2 3;
      2 3 0 1 3 2 1 2;
      4 0 0 0 4 3 2 1;
      1 2 0 5 0 1 2 3;
      0 2 0 2 10 0 1 2;
      0 2 0 2 0 5 0 1;
      6 0 5 10 0 1 10 0];
  
size=length(D(1,:));
  
for i=1:size
    for j=1:size
        if (i>j)       
            D(i,j)=D(j,i);
        end % end if
    end %end for
end % end for

end % end function dist
%******************************************************************
%******************************************************************
function shouldGo = ShouldAntGo(ant,AntSystem,n,par)
shouldGo = 1;

if all(AntSystem.Tours(ant,n,:))
    shouldGo=0;
end % end if

end %end function ShouldAntGo