function [CI,K,kp,D]=kmeans_ci_LA(X,stand,prop,nclus,nsim);
% [CI,K]=kmeans_ci_LA(X,stand,prop,nclus,nsim);
%
% This function makes a k-means clustering of the rows of the input matrix
% 'X' and computes the classifiability index for a partition of the
% rows of the input matrix 'X' using a dynamical clustering algorithm
% (k-means) into 'nclus' clusters. There are 'nsim' different partitions
% that are compared. The clustering is performed in EOF-space taking the
% leading EOF accouting for 'prop' proprtion of the total variance. The
% algorithm computes for each cluster the mean anomaly correlation between
% its centroid and the other clusters.
%
% Input
% 'X': input matrix to be classified (rows are classified)
% 'stand': option for standardising the matrix (by columns) = 'm' for just
% removing the long-term mean and ='s' for standardising to zero mean and
% unit variance. if stand is empty, the matrix is not standardized.
% 'prop': scalar giving the proportion of variance to be retained in the
% EOF pre-filtering step. if prop is empty, no EOF-prefiltering is done
% 'nclus': number of cluster
% 'nsim' : number of different partition simulated (typically nsim= 50-100)
%
% Outpurs
% 'CI': scalar giving the classifiability index
% 'K': row vector giving the cluster for each row of input matrix 'X'.
%
% ref. Michelangeli et al., JAS, 1995 (1230-1246)
%
% Vincent Moron
% June 2006
%
% L.Agel
% changed to calculate CI more consistently to Michelangeli
% documented and changed some index names
% singleton problem (single item per cluster) fixed as in kmeans_ci2.m
% random generator as for kmeans_ci2.m
% also calculates euclidean distance to centroid for each pattern
%D returns euclidiean distance for EACH simulation
%KS returns cluster assignment for each simulation

stream=RandStream('mrg32k3a');

if ~isempty(stand);
    X=stan(X,stand);
end

if ~isempty(prop);
    [U,S,V]=svd(X,0);
    s=diag(S).^2;
    sc=s/sum(s);
    a=find(cumsum(sc)> prop);
    a=a(1);
    X=U(:,1:a)*S(1:a,1:a);
end

[nt,nx]=size(X);
    
D=nan(nt,nsim);
MC=nan(nx*nclus,nsim);  %to store centroids of each partition/cluster
kp=nan(nt,nsim);        %to store results of kmeans for each partition 
for p=1:nsim; %for each partitioning
    [KK,~,~,DD]=kmeans(X,nclus,'Maxiter',1000,'EmptyAction','drop'); %'drop'==NaN or 'singleton'
    kp(:,p)=KK;
    %calculate centroid of each cluster
    for k=1:nclus;
        a=find(kp(:,p)==k);
        %added by L.Agel for case of singleton (as done for kmeans_ci2)
        if length(a)==1
            mean_cluster(k,:)=X(a,:);
        else
            mean_cluster(k,:)=nanmean(X(a,:));  %now size[nclus,nx]
        end
    end
    %make anomalies and divide by sqrt(sum(anomalies^2)) -- use this
    %instead of stan (problems with dividing by sample size)
    my=ones(nx,1)*nanmean(mean_cluster');
    mean_cluster2=mean_cluster'-my;
    sty=ones(nx,1)*sqrt(nansum(mean_cluster2.^2));
    mean_cluster2=mean_cluster2./sty;
    %lay them into MC with one column for each partition
    MC(:,p)=mean_cluster2(:);  %size(nx*numclus,nsim)
    %lay in distance matrix for D output
    for k=1:nclus
        a=find(KK==k);
        D(a,p)=DD(a,k); %DD(nt,nclus], choose the one it belongs to
    end
end

% Calculate anomaly correlation coefficients (ACCs) between two partition clusters Pi,Qj
c=nan(nsim); %to store chosen ACC to represent correlation between every partition pair P,Q
for p=1:nsim; %for each partition P
    a=find(p~=[1:nsim]); 
    for q=1:numel(a)  %compare to partition Q
        A=nan(nclus); 
        for i=1:nclus %compare cluster Pi  
            for j=1:nclus %to clusters Qj
                %calculate ACC
                 Pi=MC(((i-1)*nx)+1:i*nx,p);
           	     Qj=MC(((j-1)*nx)+1:j*nx,a(q));
                 A(i,j)=Pi'*Qj;              
            end
        end 
        %first pick best correlation of each cluster Pi to any clusters Qj in Q
        %then pick minimum of those to represent the correlation of P to Q
        c(p,a(q))=nanmin(nanmax(A')); 
    end
end

%CI = average of all partition correlations not including c(P,P)
CI=nanmean(nanmean(c));

%best partition
part=find(nanmean(c)==max(nanmean(c)));
if length(part) > 1; part=part(1); end
K=kp(:,part); %partition that maximize the ACC with the nsim-1 other partitions

% %now compute the euclidean distance of each point to the cluster mean
% cc=nan(nclus,nx);
% dd=zeros(nt,nx);
% for k=1:nclus
%     cc(k,:)=mean(X(K==k,:));
%     for j=1:nx
%         dd(K==k,j)=(X(K==k,j)-cc(k,j)).^2;
%     end
% end
% D=nansum(dd,2);