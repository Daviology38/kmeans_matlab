%do kclust analysis on amultiple variables -- 
%raw fields
%%
clear all; close all;

outdir='allthree';
if ~isdir(outdir)
    mkdir(outdir); 
end

retain=.85;
maxclust=10;
maxclust2=10;
nsim=100;  %for Monte Carlo
stream=RandStream('mrg32k3a');

%get grid
lats = ncread("CNRM-CM6-1-HR_zg.nc",'lat');
lons = ncread("CNRM-CM6-1-HR/CNRM-CM6-1-HR_zg.nc",'lon');
[LON,LAT]=meshgrid(lons,lats);

h500 = zeros(3680,24,37);
mslp = zeros(3680,24,37);
u850 = zeros(3680,24,37);
v850 = zeros(3680,24,37);

h500m = zeros(60,141);
mslpm = zeros(81,121);
u850m = zeros(81,121);
v850m = zeros(81,121);
count = 1;



z = ncread("CNRM-CM6-1-HR/CNRM-CM6-1-HR_zg.nc",'zg');
%start = 10645; %Spring start
start = 10829; %Fall start
j = start;
year = 1979;
count = 1;
counter = 1;
while j <= length(z)
       htmp = z(:,:,1,j);
        for ii = 1:60
            for jj = 1:141
                h500(counter,ii,jj) = htmp(jj,ii);
                h500m(ii,jj) = h500m(ii,jj) + htmp(jj,ii);
            end
        end
        j = j + 1;
        count = count + 1;
        counter = counter + 1;
        if(count == 92)
            %j = j+273;
            j = j+274;
            year = year + 1;
            count = 1;
            if(mod(year,4) == 0)
                j = j + 1;
            end
        end
        
end

count = count - 1;

h500m = h500m ./ 3277; %Fall season

% Remove the overall mean
for i = 1:length(h500)
   h500(i,:,:) = squeeze(h500(i,:,:)) - h500m(:,:);
   mslp(i,:,:) = squeeze(mslp(i,:,:)) - mslpm(:,:);
   u850(i,:,:) = squeeze(u850(i,:,:)) - u850m(:,:);
   v850(i,:,:) = squeeze(v850(i,:,:)) - v850m(:,:);
end

%area weight
for l=1:length(h500)
    h500aw(l,:,:)=squeeze(h500(l,:,:)).*sqrt(cos(pi*LAT/180.));
    mslpaw(l,:,:)=squeeze(mslp(l,:,:)).*sqrt(cos(pi*LAT/180.));
    u850aw(l,:,:)=squeeze(u850(l,:,:)).*sqrt(cos(pi*LAT/180.));
    v850aw(l,:,:)=squeeze(v850(l,:,:)).*sqrt(cos(pi*LAT/180.));
end


%arrange as time x space
F=ones(size(LON));
u850u=eof_map2mat(F,u850aw);
v850u=eof_map2mat(F,v850aw); 
h500u=eof_map2mat(F,h500aw); 
mslpu=eof_map2mat(F,mslpaw); 

%put winds together
multiu= [h500u mslpu u850u v850u];

%standardize values at each grid point
multius=stan(multiu,'s');

multius(isnan(multius))=0;
multius(multius==-32767)=0;

%eof going in to reduce data size
[U,S,V]=svd(multius,0);
s=diag(S).^2;
sc=s/sum(s);
a=find(cumsum(sc)>retain);
a=a(1);
tmpu=U(:,1:a)*S(1:a,1:a);
nr=size(tmpu,1); nc=size(tmpu,2);
%    
% %sanity check
res=input(sprintf('Perform k-means on dates,  variables? (Y/N):'),'s');
if ~strcmp(upper(res),'Y')
    display('...ending')
    return;
end
%%
%determine optimum clusters from Classifiability Index (CI)
if ~exist(sprintf('%s/CI_results.mat',outdir),'file') 
    for i=1:maxclust
        display(sprintf('k=%d',i));
        [CI(i),K(:,i),D(:,:,i)]=kmeans_ci2(tmpu,[],[],[],i,100);
    end
    save(sprintf('%s/CI_results',outdir),'CI','K','retain','D');
else
    tmp=load(sprintf('%s/CI_results',outdir));
    CI=tmp.CI;
    K=tmp.K; 
    D=tmp.D;
end

figure
hold on
plot(2:maxclust,CI(2:maxclust));
hold on
ylabel('CI');
xlabel('Cluster');
title('CI values, multi');
print(gcf,'-dpng',sprintf('%s/plot_ci.png',outdir));

%%
%create nsim random noise samples (only do this once - use for all cluster sizes)
%from the standardized anomalies
if ~exist(sprintf('%s/random.mat',outdir),'file')
    display('creating random red noise series');
    rsims=nan(nr,nc,nsim); %nc is reduced through EOF
    for i=1:nc
        p=ar1rand(tmpu(:,i),nsim);
        %p=ebisuzaki(tmpu(:,i),nsim,-1);
        rsims(:,i,:)=p;
    end
    save(sprintf('%s/random',outdir),'rsims');
else
    tmp=load(sprintf('%s/random',outdir));
    rsims=tmp.rsims;
end

%for each cluster size, do nsim times using noise
CIcis=nan(maxclust2,nsim);
for i=1:maxclust2
    if ~exist(sprintf('%s/CICI_%d.mat',outdir,i),'file')
        for j=1:nsim
            disp(sprintf('Cluster %d Simulation %d',i,j));
            sim=squeeze(rsims(:,:,j));
            CIci(1,j)=kmeans_ci2(sim,[],[],[],i,100);
        end
        save(sprintf('%s/CIci_%d',outdir,i),'CIci');
    end
    tmp=load(sprintf('%s/CIci_%d',outdir,i));
    CIcis(i,:)=tmp.CIci;
    cisort=sort(CIcis(i,:));
    citop(i)=cisort(.90*nsim);  %one-sided 90% confidence interval
    cibottom(i)=cisort(1);
    plot([i i],[cibottom(i) citop(i)],'-r');
    drawnow
end

%formal picture
figure
fill([1:maxclust2 maxclust2:-1:1],[cibottom(1:maxclust2) citop(maxclust2:-1:1)],rgb('LightGray'),'Edgecolor','none');
hold on
plot(1:maxclust,CI(1:maxclust),'LineWidth',4);
ylabel('CI');
xlabel('Cluster');
title('CI values with confidence intervals, multi');
print(gcf,'-dpng',sprintf('%s/plot_cici.png',outdir));

%do dissimilarity index
for k=1:maxclust
    for i=1:k
        f=K(:,k)==i;
        tmpxmu(i,:)=nanmean(tmpu(f,:,:));
    end
    ACC=nan(k,k);
    for i=1:k
        for j=1:k
            Pi=tmpxmu(i,:);
            Pj=tmpxmu(j,:);
            ACC(i,j)=corr(Pi',Pj');
        end
    end
    minACC(k)=nanmin(nanmin(ACC.^2));
end
figure
plot(1-minACC);
xlabel('Clusters');
ylabel('1-ACC^2');
title('Dissimilarity index multi');
print(gcf,'-dpng',sprintf('%s/kclust_dissim.png',outdir));

% %try Riddle et al. 2012 technique
%fourth way, figure out distance per map, which is sum of
%values in D, then take 90% of those for each cluster
DD=squeeze(sum(D,2));  %now DD(nr,k)
R=nan(maxclust,maxclust);
for k=1:maxclust
    for j=1:k
        f=K(:,k)==j;
        krads=D(f,k); 
        krads=sort(krads(:));
        R(k,j)=krads(floor(.9*numel(krads)));
    end
end
for k=1:maxclust
    sigma(k)=nansum((R(k,:)/R(1,1)).^nc);
end
figure
plot(sigma)
title('Volume ratio index, 500 hPa geop');
print(gcf,'-dpng',sprintf('%s/volume_ratio_index.png',outdir));

%%
%for response
figure
subplot_tight(3,3,1,[.09 .09]);
hold on
fill([1:maxclust2 maxclust2:-1:1],[cibottom(1:maxclust2) citop(maxclust2:-1:1)],rgb('LightGray'),'Edgecolor','none');
plot(CI(1:maxclust));
%ylabel('CI','fontsize',14);
xlabel('Cluster','fontsize',12);
ylim([0.5 1]);
xlim([1 10]);
set(gca,'Xtick',[1:9]);
set(gca,'Ytick',[0.5:.1:1]);
title('Classifiability Index','fontsize',12);
%
subplot_tight(3,3,2,[.08 .08]);
plot(1-minACC);
xlabel('Cluster','fontsize',12);
%ylabel('DI','fontsize',14);
xlim([1 10]);
set(gca,'Xtick',[1:9]);
title('Dissimilarity Index','fontsize',12);
%
subplot_tight(3,3,3,[.08 .08]);
plot(sigma)
xlabel('Cluster','fontsize',12);
%ylabel('VRI','fontsize',14);
ylim([0 1]);
xlim([1 10]);
set(gca,'Xtick',[1:9]);title('Volume Ratio Index','fontsize',12);
print(gcf,'-dpng',sprintf('%s/indices_multi.png',outdir));

for i = 1:10
figure
[silh3,h] = silhouette(tmpu,K(:,i));
mean(silh3>0)
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value'
ylabel 'Cluster'
end