function X=nan_std(Y,missing);

% X=nan_std(Y,missing);
%
% This function replaces nanstd when it is not available.
%
% Vincent Moron
% July 2006

if nargin==1;
    missing=NaN;
end

if missing ~= NaN;
    Y(find(Y==missing))=NaN*ones(size(find(Y==missing)));
end

[nr,nc]=size(Y);
nnans=~isnan(Y);
Y2=Y;
Y2(~nnans)=0;
n=sum(nnans);
n(find(n==0))=NaN;
YM=ones(nr,1)*nan_mean(Y);
X=sqrt(sum(((Y2-YM).^2).*nnans)./(n-1));

