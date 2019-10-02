%%%%%%%%%% Option parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = 50; % Value of the underlying
K = 50; % Strike (exercise price)
r = 0.05; % Risk free interest rate
sigma = 0.2; % Volatility
T = 3; % Time to expiry
%%%%%%%%%% Method Paramaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 80; % Number of asset mesh points
N = 200; % Number of time mesh points
Szero = 0; % Specify extremes
Smax= 150;
%%%%%%%%%% Setup Mesh and BCs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solution_mesh=zeros(N+1,M+1); % Create Mesh
Smesh=0:(Smax/M):Smax; % Mesh of equally spaced S values
Tmesh=T:-(T/N):0; % Mesh of equally spaced T values
dt=T/N; % Specify timestep, dt
solution_mesh(1,:)= max(K-Smesh,0); % Payoff BC (Put)
solution_mesh(:,1)= K*exp(-r*(T-Tmesh)); % BC at S=0
solution_mesh(:,M+1)= 0; % BC at S=M
A = @(i) 0.5*dt*(r*i-sigma^2*i^2);
B = @(i) 1 + (sigma^2*i^2 + r)*dt; % Define the functions A, B & C
C = @(i) -0.5*dt*(sigma^2*i^2+r*i);
%%%%%%%%%% Construct Tridiagonal Matrix, Tri %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Acoeffs = zeros(M+1,1);Bcoeffs = zeros(M+1,1);Ccoeffs = zeros(M+1,1);
for i=1:M+1
Acoeffs(i) = A(i-1); Bcoeffs(i) = B(i-1); Ccoeffs(i) = C(i-1);
end
Tri=diag(Acoeffs(2:end),-1)+diag(Bcoeffs)+diag(Ccoeffs(1:end-1),+1);
Tri_Inv=inv(Tri); % Compute inverse
%%%%%%%%%% Implicit Euler Iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:N
temp=zeros(M+1,1);
temp(1)=A(0)*solution_mesh(j+1,1);
temp(end)=C(M)*solution_mesh(j+1,M+1); % Boundary terms
RHS=solution_mesh(j,:)'-temp;
temp=Tri_Inv*RHS;
solution_mesh(j+1,(2:end-1))=temp(2:end-1);
end

%Save file to CSV
%filename = 'FuzzySets.csv'; 
%dlmwrite(filename, solution_mesh,'-append', 'delimiter', ',');

%filename = 'CrispInputs.csv'; 
%dlmwrite(filename, Smesh,'-append', 'delimiter', ',');

%filename = 'FuzzyRules.csv'; 
%dlmwrite(filename, Tmesh,'-append', 'delimiter', ',');

%Apply fuzzyLogic
%FuzzySet = FLP_LoadFuzzySets('FuzzySets.csv');
%CrispInput = FLP_LoadCrispInput(FuzzySet,'CrispInputs.csv');
%FuzzyRules = FLP_LoadFuzzyRules(FuzzySet,'FuzzyRules.csv');

%AntMemberGrades = FLP_Fuzzification(FuzzySet, CrispInput);
%ConsqMemberGrades = FLP_FuzzyRuleEval(AntMemberGrades,FuzzyRules);
%[CrispOutput, OutputMF, X] = FLP_DeFuzzification(ConsqMemberGrades, FuzzySet, 100);
%plotdata = FLP_3DMFplot( FuzzySet, FuzzyRules, 1, 2, 20 );

load('FuzzySets.csv');
load('FuzzyRules.csv');
FIS_type = 'mamdani'; % You can also set it to 'sugeno'

% Ignore the last column if it consists of only zeros
if sum(FuzzySets(:,end)) == 0
    Xin = FuzzySets(:,1:end-1); % Input
else
    Xin = FuzzySets; % Input
end

l = size(FuzzyRules,2);
noOfPoints = size(FuzzySets,1);
% Obtain the outputs corresponding to the inputs
if l == noOfPoints
    Xout = FuzzyRules;
else if (numel(FuzzyRules) == size(FuzzySets,1))
    Xout = repmat(FuzzyRules(1,:),1,3);
    end
end

% fcm options
n_clusters = 5; % No of Clusters to be used in the fcm
fcm_opt = NaN(4,1);
fcm_opt(4) = 0;

% Generate the Fuzzy Inference System
FIS_mat = genfis3(Xin, Xout, FIS_type, n_clusters, fcm_opt);

% Plot membership functions
ROWS_PER_FIG = 2; % No. of rows in a figure
COLUMNS_PER_FIG = 3; % No. of columns in a figure
PLOTS_PER_FIG = ROWS_PER_FIG * COLUMNS_PER_FIG;
fig_pos = 1;
fig = figure;
% Maximize the figure screen
set (fig, 'Units', 'normalized', 'Position', [0,0,1,1]);
for input_no = 1:size(Xin,2)
    [x,mf] = plotmf(FIS_mat,'input',input_no);
    subplot(ROWS_PER_FIG, COLUMNS_PER_FIG, fig_pos), plot(x,mf)
    xlabel(sprintf('Membership Functions for Input %d', input_no))
    
    if mod(fig_pos,PLOTS_PER_FIG) == 0
        fig_pos = 1;
        fig = figure;
        set (fig, 'Units', 'normalized', 'Position', [0,0,1,1]);
    else
        fig_pos = fig_pos + 1;
    end
end

% Surface plot (Change the second parameter to try out different input
% combinations. Example: If you want to see the plots of inputs 5 and 47,
% set the second parameter as [5,47].
figure;
gensurf(FIS_mat, [3,80],1);
% Uncomment the next line to open surface editor
% surfview(FIS_mat)

% Rule Viewer
ruleview(FIS_mat);

%Find 5 clusters using fuzzy c-means clustering.
figure;
[centers,U,obj_fcn] = fcm(FuzzySets,n_clusters);

%Classify each data point into the cluster with the largest membership
%value.
maxU = max(U);
index1 = find(U(1,:) == maxU);
index2 = find(U(2,:) == maxU);
index3 = find(U(3,:) == maxU);
index4 = find(U(4,:) == maxU);
index5 = find(U(5,:) == maxU);

%Plot the clustered data and cluster centers.
plot(FuzzySets(index1,1),FuzzySets(index1,5),'ob')
hold on
plot(FuzzySets(index2,1),FuzzySets(index2,5),'or')
hold on
plot(FuzzySets(index3,1),FuzzySets(index3,5),'og')
hold on
plot(FuzzySets(index4,1),FuzzySets(index4,5),'oy')
hold on
plot(FuzzySets(index5,1),FuzzySets(index5,5),'ob')
hold on
plot(centers(1,1),centers(1,5),'xb','MarkerSize',15,'LineWidth',3)
plot(centers(2,1),centers(2,5),'xr','MarkerSize',15,'LineWidth',3)
plot(centers(3,1),centers(3,5),'xg','MarkerSize',15,'LineWidth',3)
plot(centers(4,1),centers(4,5),'xy','MarkerSize',15,'LineWidth',3)
plot(centers(5,1),centers(5,5),'xb','MarkerSize',15,'LineWidth',3)
hold off

save solution_mesh.dat solution_mesh -ascii
fcmdemo
mesh(Smesh,Tmesh,solution_mesh)
xlabel('S');ylabel('t');
zlabel('V(S,t)')
%%%%%%%%% Extract Desired Values Using Interpolation %%%%%%%%%%%%%%%%%%%%%%
interp1(Smesh,solution_mesh(N+1,:),S)
