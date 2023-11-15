function [Z]=scale_mean_var(X,Y);

% [Z]=scale_mean_var(X,Y);
%
% This function scales the columns of the input matrix 'X', so that its 
% mean and variance match the ones of the input matrix 'Y'.
%
% Input
% 'X' : matrix of real number to be scaled
% 'Y' : matrix of real number used in the scaling
%
% Output
% 'Z' : matrix of real number (the mean and variance of the columns of 'Z'
% are equal to the ones of 'Y'
%
% Vincent Moron
% Nov 2005

[nr,nc]=size(X);
my=mean(Y);
sy=std(Y);
mx=mean(X);
sx=std(X);

Z=(X*diag(sy./sx));
dm=mean(Z)-my;
dm=ones(nr,1)*dm;
Z=Z-dm;