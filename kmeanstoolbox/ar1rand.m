function Y=ar1rand(X,nsim);

% Y=ar1rand(X,nsim);
%
% Red-noise random simulation of 'nsim' time series of the
% same length as the input vector 'X' and having the same
% one-order auto-correlation and mean and variance as 'X'
%
% Vincent Moron
% June 2006

X=X(:);
n=length(X);
c=corrcoef(X(1:n-1),X(2:n));
c=c(1,2);
d=1-c;
Y=randn(1,nsim);
Z=randn(n,nsim);
Z=scale_mean_var(Z,copy(X(2:n)',nsim)');
for j=2:n;
    Y(j,:)=(Y(j-1,:).*c)+(Z(j,:).*d);
end
Y=scale_mean_var(Y,copy(X',nsim)');
