function [X]=copy(Y,N);

% [X]=copy(Y,N);
%
% This function copies 'N' times a matrix 'Y'.
%
% Input
% 'X' : matrix of real number
% 'N' : scalar integer giving the number of copies
%
% Output
% 'Y' : matrix 'X' copied 'N' times (the copy is done in row)
%
% Vincent Moron
% July 1997

[NR,NC]=size(Y);
X=Y';
X=X(:);
X=X*ones(1,N);
X=reshape(X,NC,NR*N);
X=X';
