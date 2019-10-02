%%%%%%%%%% Option parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = 3; % Value of the underlying
K = 1; % Strike (exercise price)
r = 0.05; % Risk free interest rate
sigma = 0.25; % Volatility
T = 3; % Time to expiry
%%%%%%%%%% Binomial Method Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 50; % Number of time steps (length of binomial tree)
dt = T/N; % Length of each time-step (automatically computed)
A = 0.5*(exp(-r* dt)+exp((r+sigma^2)*dt));
u = A + sqrt(A^2-1); d = A - sqrt(A^2-1);
p = (exp(r*dt)-d)/(u-d);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Underlying values at time T
P=pascal(N+1); % Pascal’s triangle
binomial_coeffs=zeros(N+1,1); % Preallocate vector for Binomial coeffs
% Extract Binomial Coefficients
for i=1:N+1
binomial_coeffs(i)=P(N+2-i,i);
end
%binomial_coeffs % This is for debugging, can comment out
%max(binomial_coeffs)
% Construct matrix for the cumulative paths (which we will later take the
% average of)
% First work out the width of the matrix (cols):
if mod(N,2)==0
cols=N/2+1;
else
cols=(N-1)/2+1;
end
% We can now preallocate the arrays:
paths=zeros(N+1,cols);
max_vals=zeros(N+1,cols)
% Now fill this matrix with the (truncated) binomial coeffients:
% Also create a matrix which holds the maximum values of the paths these

% correspond to the coefficients which were created in in ’paths’:
for j=1:cols
for i=j:N+2-j
if j==1
paths(i,j)=1;
max_vals(i,j)=u^(N+1-i);
else
paths(i,j)=binomial_coeffs(j)-sum(paths(i,1:j));
max_vals(i,j)=u^(max_vals(i,j-1)-1);
end
end
end
%paths
%max_vals
option_vals=max(max_vals*S-K,0); % Option values for each maxima
% final_vals is the (.*) product of option_vals and paths, summed across
% the columns and then (mean) averaged by divig through with binomial_coeffs:
Tree=sum(option_vals.*paths,2)./binomial_coeffs;
Back = p*eye(N+1,N+1) + (1-p)*diag(ones(N,1),1);
Back = sparse(Back); % Define matrix to track back through tree
%%%%%%%%%% Track back through the tree to time t=0 %%%%%%%%%%%%%%%%%%%%%%%%
for i = N:-1:1
Tree = Back(1:i,1:i+1)*Tree;
end
%%%%%%%%%% Discount under risk-neutral assumption %%%%%%%%%%%%%%%%%%%%%%%%%
Lookback_Value = exp(-r*T)*Tree;
display(Lookback_Value)

