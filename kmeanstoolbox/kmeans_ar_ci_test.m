function [CI]=kmeans_ar_ci_test(X,stand,prop,nclus,nsim);

% [CI]=kmeans_ar_ci_test(X,stand,prop,nclus,nsim);
%
% This function computes a red-noise test for testing the significance of 
% the classifiability index used to determine the best number in the
% dynamical clustering scheme.
%
% Inputs
% 'X': input matrix to be classified (rows are classified)
% 'stand': option for standardising the matrix (by columns) = 'm' for just
% removing the long-term mean and ='s' for standardising to zero mean and
% unit variance. if stand is empty, the matrix is not standardized.
% 'prop': scalar giving the proportion of variance to be retained in the
% EOF pre-filtering step. if prop is empty, no EOF-prefiltering is done
% 'nclus': number of cluster
% 'nsim' : number of different partition simulated and different
% simulations (typically nsim= 50-100)
%
% Outpurs
% 'CI': scalar giving the classifiability index for the red-noise.
%
% ref. Michelangeli et al., JAS, 1995 (1230-1246)
%
% Vincent Moron
% June 2006
%
% L. Agel Jan 2015
% check for empty EOF filtering

[r,c]=size(X);

if ~isempty(stand)
    X=stan(X,stand);
end

%L.Agel added check for empty prop, allow for full analysis, not reduced
if ~isempty(prop)
    [U,S,V]=svd(X,0);
    s=diag(S).^2;
    sc=s./sum(s);
    a=find(cumsum(sc) > prop);
    a=a(1);
    PC=U(:,1:a)*S(1:a,1:a);
    N=a;
else
    PC=X;
    N=size(X,2);
end

for i=1:N;
    XR(:,(i:N:N*nsim))=ar1rand(PC(:,i),nsim); 
end

for i=1:nsim;
    display(['simulation # ',num2str(i)]);
    [CI(i),k]=kmeans_ci(XR(:,((i-1)*N)+1:i*N),[],[],nclus,nsim);
end
