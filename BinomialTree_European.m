%%%%%%%%%% Option parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = 55; % Value of the underlying
K = 50; % Strike (exercise price)
r = 0.05; % Risk free interest rate
sigma = 0.25; % Volatility
T = 3; % Time to expiry
%%%%%%%%%% Binomial method parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 256;
dt = T/N; A = 0.5*(exp(-r* dt)+exp((r+sigma^2)*dt));
d = A - sqrt(A^2-1); u = A + sqrt(A^2-1);
p = (exp(r*dt)-d)/(u-d);
%%%%%%%%%% Construct the tree %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tree=S; % Initialise tree
for i=1:N
mult=u*ones(i+1,1); % Create array of ’up multiplications’
mult(end)=d; % Make the last entry a ’down multiplication’
Tree(i+1)=Tree(i);
if i==1
Tree=Tree'; % Make sure T is in column form
end
Tree=Tree.*mult;
end
Tree=max(K-Tree,0) % Compute Put option values at t=T
Back = p*eye(N+1,N+1) + (1-p)*diag(ones(N,1),1);
Back = sparse(Back); % Define matrix to track back through tree

%%%%%%%%%% Track back through the tree to time t=0 %%%%%%%%%%%%%%%%%%%%%%%%
for i = N:-1:1
Tree = Back(1:i,1:i+1)*Tree;
end

%%%%%%%%%% Discount under risk-neutral assumption %%%%%%%%%%%%%%%%%%%%%%%%%
Option_Value = exp(-r*T)*Tree;
display(Option_Value)
save BinomialTreeEuropean.dat Tree -ascii
fcmdemo
