%%%%%%%%%% Option parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = 50; % Value of the underlying
K = 50; % Strike (exercise price)
r = 0.05; % Risk free interest rate
sigma = 0.25; % Volatility
T = 3;
%%%%%%%%%% Binomial method parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 512;
dt = T/N; A = 0.5*(exp(-r* dt)+exp((r+sigma^2)*dt));
up = A + sqrt(A^2-1); down = A - sqrt(A^2-1);
p = (exp(r*dt)-down)/(up-down);
%%%%%%%%%% Construct the tree %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tree = S*up.^((N:-1:0)').*down.^((0:N)');
% Compute Put option values at t=T
Tree=max(K-Tree,0);
% Matrix to compute expectations
Back = p*eye(N+1,N+1) + (1-p)*diag(ones(N,1),1);
%%%%%%%%%% Track back through the tree to time t=0 %%%%%%%%%%%%%%%%%%%%%%%%
for i = N:-1:1
discounted = Back(1:i,1:i+1)*Tree*exp(-r*dt); % One step back, discounted
sharevals = S*up.^((i-1:-1:0)').*down.^((0:i-1)'); % Share values at time step
Tree=max(discounted,K-sharevals); % American Exercise Comparison
end
%%%%%%%%%% Print Option Value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Option_Value = Tree;
display(Option_Value)
