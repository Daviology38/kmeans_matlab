function [X]=ebisuzaki(x,nsim,value);

% [X]=ebisuzaki(x,nsim,value);
%
% This function creates 'nsim' random time series that have the same power
% spectrum as the original time series 'x' but random phases. The power
% spectrum of 'x' and 'X' is the same (see abs(fft(x)) and abs(fft(X));
% The 'X' time series has also the same variance as 'x'. (see Ebisuzaki W. (1997) ;
% A method to estimate the statistical significance of a correlation when 
% the data are serially correlated. J Climate, 10, 2147-2153).
%
% Input
% 'x' : vector of real number containing the time series to match
% 'nsim' : integer number giving the number of simulation
% 'value' : real number to initiate the random sequence (if > 0, the seeds
% is initiated to the value; otherwise, it is changed from the clock)
%
% Output
% 'X' : matrix of real number (the same number of rows as length of 'x' and
% 'nsim' columns.
%
% Vincent MORON
% Mars 2002

if(value>0);
    rand('seed',value);
    randn('seed',value);
else
%    rand('seed',sum(100*clock));
%    randn('seed',sum(100*clock));
   stream=RandStream('mrg32k3a');
end

n=length(x);
n2=floor(n/2);
x=x(:);
x=stan(x,'m');
y=fft(x); % DFT de la s�rie originale
mod=abs(y); % puissance spectrale observ�e

for i=1:nsim;
   if n/2==n2 % cas des s�ries paires
      an1=angle(y(2:n2));
      tt=randn(n2-1,1);
      anr1=(tt.*2.*pi)-pi;
      anr1=[0;anr1(:);0;flipud(anr1(:).*(-1))];
      recf=mod.*exp(sqrt(-1)*anr1); % DFT des s�ries permut�es 
      X(:,i)=real(ifft(recf));
   else % cas des s�ries impaires
      an1=angle(y(2:n2+1));
      tt=randn(n2,1);
      anr1=(tt.*2.*pi)-pi;
      anr1=[0;anr1(:);flipud(anr1(:).*(-1))];
      recf=mod.*exp(sqrt(-1)*anr1); % DFT des s�ries permut�es
      X(:,i)=real(ifft(recf));
   end
end

% standardisation des s�ries permut�es relativement
% � la variance de la s�rie observ�e

X=X*diag(std(x)./std(X));