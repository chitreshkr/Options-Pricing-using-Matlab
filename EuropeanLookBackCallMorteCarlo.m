%%%%%%%%%% Option parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = 3; % Value of the underlying
K = 1; % Strike (exercise price)
r = 0.05; % Risk free interest rate
sigma = 0.25; % Volatility
T = 3; % Time to expiry
%%%%%%%%%% Monte-Carlo Method Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% randn(’state’,0)
M=5000; % Number of Monte-Carlo trials
n=1000; % Set number of steps to compute at in [0,T]
dt=T/n; % Compute the time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dW=sqrt(dt)*randn(M,n); % Generate array of brownian movements drawn from N(0,dt^2)
W=cumsum(dW,2); % Sum the array cummulativley to obtain Wiener process
t=0:dt:T; % Array of equal time steps
W=[zeros(M,1),W]; % Set first values to be zero
tt=repmat(t,M,1); % Create matrix for t vals
% Compute the asset path using the solution of the asset path SDE
asset_path=S*exp((r-0.5*sigma^2)*tt+sigma*W);
% Compute the maximum value the asset reaches over the life of the option:
max_vals=max(asset_path,[],2);
% Evaluate the fixed strike lookback call option in each case:
option_values=max(max_vals-K,0);
% Discount under risk-neutral assumption
present_vals=exp(-r*T)*option_values;
call_value=mean(present_vals);
display(call_value)
