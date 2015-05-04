function y = mvnrnd_spinlab(mu,Sigma,L)
% Radu David April 2015
% Make multivariate random numbers given mu and Sigma
% If Sigma is singular, a small diagonal offset is added to do Cholesky
% factorization
% Need to pass in a vector mu and a covariance matrix Sigma with dimensions
% of square


narginchk(2,3); % check number of arguments

if nargin==2 %if only two arguments given, generate one random evaluation
    L = 1;
end


n = length(mu);
% generate uniform values for x
x = zeros(n,L);
%using Cholesky decomposition
n = size(Sigma,1);
for i = 1:n
    x(i,:) = randn(1,L);
end
% check if Sigma is singular, if yes, add a small diagonal offset to be
% able to run Cholesky factorization
tol = 10*eps(max(abs(diag(Sigma))));
if (all(eig(Sigma))<tol)
    Sigma = Sigma + 1e-40*eye(n);
end
ch = chol(Sigma);
y = ch'*x+mu'*randn(1,L);
y = y';
