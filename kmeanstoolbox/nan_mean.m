function X=nan_mean(Y,missing);

% X=nan_mean(Y,missing);
%
% This function replaces nanmean when it is not available.
%
% Vincent Moron
% July 2006

if nargin==1;
    missing=NaN;
end

if missing ~= NaN;
    Y(find(Y==missing))=NaN*ones(size(find(Y==missing)));
end

nans=isnan(Y);
Y(nans)=0;
n=sum(~nans);
n(find(n==0))=NaN;
X=sum(Y)./n;
