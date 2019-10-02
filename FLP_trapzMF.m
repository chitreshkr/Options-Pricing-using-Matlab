function [ mf ] = FLP_trapzMF( x,tparms )
% FLP_trapzMF Calculates the value of the trapezoidal membership function at
% each point x. 
%
% Input
% x - a column vector of values
% tparms - a vector with the a, b, c, & d trapezoidal parameters
%
% Output
% a columnar vector with membership function values
%
% Author: Jim Kunce (jdk_acct@yahoo.com)

mf = zeros(size(x,1),1); % preformat the output with zeros
a = tparms(1); % get the trapezoidal parameters
b = tparms(2);
c = tparms(3);
d = tparms(4);

% Use logical indexing to identify where each x value falls within the
% trapezoidal function
idx1 = x(:,1) > a & x(:,1) < b;
idx2 = x(:,1) >= b & x(:,1) <= c;
idx3 = x(:,1) > c & x(:,1) < d;

% Replicate the scalar parameters to vectors to facilitate a "vectorized"
% calculation of the membership function
a = repmat(a,size(x,1),1);
b = repmat(b,size(x,1),1);
c = repmat(c,size(x,1),1);
d = repmat(d,size(x,1),1);

% Calculate and assign the membership function for each range
res = (x-a) ./ (b-a); mf(idx1,1) = res(idx1,1);
mf(idx2,1) = 1;
res = (d-x) ./ (d-c); mf(idx3,1) = res(idx3,1);

end

