function [done]=kmeans_ar_ci_test_LA(X,stand,prop,maxclusts,nsim,serial);

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

% L.Agel Jan 2015 
% nclus changed to maxclusts (do all here within this routine)
% changed to call kmeans_ci_LA (does CI calculation more consistent with Michelangeli)
% different radomization for serial correlation
% returns CIci[nclus,nsim]

%can standardize or remove global mean or nothing

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
else
    PC=X;
end

[nr,nc]=size(PC);

%L.Agel create red noise samples here, just do once or all cluster sizes
if ~exist('random.mat','file')
    display('creating random red noise series');
    rsims=NaN*ones(nr,nc,nsim);
    for i=1:nc
        if serial==1
            p=ebisuzaki(PC(:,i),nsim,-1); 
        else
            p=ar1rand(PC(:,i),nsim);
        end
        rsims(:,i,:)=p;
    end
    save('random','rsims');
else
    tmp=load('random');
    rsims=tmp.rsims;
end

%for each cluster size, do nsim times using noise
for i=1:maxclusts
    if ~exist(sprintf('CICI_%d.mat',i),'file')
        for j=1:nsim
            disp(sprintf('Cluster %d Simulation %d',i,j));
            sim=squeeze(rsims(:,:,j));
            cis(1,j)=kmeans_ci_LA(sim,[],[],i,100);
        end
        save(sprintf('CIci_%d',i),'cis');
    end
end

done=1;
