%%%%%%%%%% Option parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = 55; % Value of the underlying
K = 50; % Strike (exercise price)
r = 0.05; % Risk free interest rate
sigma = 0.25; % Volatility
T = 3; % Time to expiry
%%%%%%%%%% Monte-Carlo Method Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% randn(’state’,0) % Repeatable trials on/off
M=1e7; % Number of Monte-Carlo trials
%%%%%%%%%% Use final values to compute %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
final_vals=S*exp((r-0.5*sigma^2)*T + sigma*sqrt(T)*randn(M,1));
option_values=max(K-final_vals,0); % Evaluate the Put option options
present_vals=exp(-r*T)*option_values; % Discount under r-n assumption
int=1.96*std(present_vals)/sqrt(M); % Compute confidence intervals
put_value=mean(present_vals); % Take the average
display(put_value)
display([put_value-int put_value+int])