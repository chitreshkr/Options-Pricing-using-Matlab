%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T _ _ _ _ _ _ _ _ _ _ _ _ _ _
% V(i,j) -> V(asset,time) | |
% t | |
% i = 0,1,2....M bcs at i=0,M | |
% j = N,N-1,...0 bc at j=N |___________________________|
% Szero S Smax
%%%%%%%%%% Option parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = 50; % Value of the underlying
K = 50; % Strike (exercise price)
r = 0.05; % Risk free interest rate
sigma = 0.25; % Volatility
T = 3;
%%%%%%%%%% Method Paramaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 100; % Number of asset mesh points
N = 400; % Number of time mesh points
Szero = 0; % Specify extremes
Smax= 150;
omega=1.2;tol=0.001; % SOR Parameters
%%%%%%%%%% Setup Mesh and BCs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solution_mesh=zeros(N+1,M+1); % Create Mesh
Smesh=0:(Smax/M):Smax; % Mesh of equally spaced S values
Tmesh=T:-(T/N):0; % Mesh of equally spaced T values
dt=T/N; % Specify timestep, dt
dS=Smax/M; % Specify asset step, dS
solution_mesh(1,:)= max(K-Smesh,0); % Payoff BC (Put)
solution_mesh(:,1)= K*exp(-r*(T-Tmesh)); % BC at S=0
solution_mesh(:,M+1)= 0; % BC at S=M
% (M+1 due to MATLAB indexes)
a = @(i) 0.5*dt*(r*i-sigma^2*i^2);
b = @(i) 1 + (sigma^2*i^2 + r)*dt; % Define the functions a, b & c
c = @(i) -0.5*dt*(sigma^2*i^2+r*i);
%%%%%%%%%% Construct Tridiagonal Matrix, Tri %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
acoeffs = zeros(M+1,1);bcoeffs = zeros(M+1,1);ccoeffs = zeros(M+1,1);
for i=1:M+1
acoeffs(i) = a(i-1);
bcoeffs(i) = b(i-1);
ccoeffs(i) = c(i-1);
end

Tri=diag(acoeffs(3:M),-1)+diag(bcoeffs(2:M))+diag(ccoeffs(2:M-1),+1);
Tri_Inv=inv(Tri); % Compute inverse
%%%%%%%%%% Implicit Euler + SOR Iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=1:N
temp=zeros(M-1,1);
temp(1)=a(0)*solution_mesh(p+1,1);
temp(end)=c(M)*solution_mesh(p+1,M+1); % Boundary terms
RHS=solution_mesh(p,2:M)'-temp;
% Gauss-Siedel SOR method: Solving Ax=b iteratively
% Note: A is (M-1)x(M-1), b is (M-1)
A=Tri;b=RHS; % Define matrix, A and RHS vector, b
x=solution_mesh(p,2:M)'; % Define x
xold=10000*x; % Initialise xold to enter loop
n=length(x);
while norm(xold-x)>tol
xold=x; % redefine xold
for i=1:n % for: rows of the matrix
if i==1
z=(b(i)-A(i,i+1)*x(i+1))/A(i,i);
x(i) = max(omega*z + (1-omega)*xold(i),K-i*dS);
elseif i==n
z=(b(i)-A(i,i-1)*x(i-1))/A(i,i);
x(i) = max(omega*z + (1-omega)*xold(i),K-i*dS);
else
z=(b(i)-A(i,i-1)*x(i-1)-A(i,i+1)*x(i+1))/ A(i,i);
x(i) = max(omega*z + (1-omega)*xold(i),K-i*dS);
end
end
end
solution_mesh(p+1,(2:end-1))=x;
end
%mesh(Smesh,Tmesh,solution_mesh)
%xlabel(’S’);ylabel(’t’);zlabel(’V(S,t)’)
%%%%%%%%% Extract Desired Value Using Interpolation %%%%%%%%%%%%%%%%%%%%%%%
interp1(Smesh,solution_mesh(N+1,:),S)

