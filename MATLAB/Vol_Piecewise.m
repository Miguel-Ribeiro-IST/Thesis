clear;
S0=167.65;
T = 28/252;
r = 0.06;
sigma0 = 0.21;
L = 56; % number of time intervals
D=252; %number of (trading) days in a year
dt = T/L;
M = 1000; % number of asset paths

A = importdata('AAPL.txt','\t',1);
B=A.data(:,:);
B=B(B(:,1)==6,:);


a=0.2154;b=   20.0000;c=  184.1392;
%a=0.1897;
%b=2.0986;
%c=179.6484;

Y = zeros(M,1);
S = S0*ones(M,L+1); % asset paths
sigma=a*ones(M,1);

for k = 2:L+1
%S(:,k)=S(:,k-1).*exp((r-0.5*sigma0^2)*dt+sigma0*sqrt(dt)*randn(M,1));
S(:,k)=S(:,k-1)+S(:,k-1)*r*dt+sqrt(dt)*sigma(:).*S(:,k-1).*randn(M,1);
if dt*D*(k-2)<=100
   sigma=arrayfun(@(x) a+b*x^2,(S(:,k)-c)/c);
elseif dt*D*k<=28
   sigma=arrayfun(@(x) a+b*x^2,(S(:,k)-c)/c);
elseif dt*D*(k-2)>150
   sigma=0*ones(M,1);
end
end

Error=0;
% Find payoff Y at expiry.
%%{
for j=1:size(B,1)
for i=1:M
Y(i,j) = max(S(i,12)-B(j,2),0);
end
end
abs(sum(exp(-r*6/D)*mean(Y)-B(:,3)'))
%}



%{
d1 = (log(S0/K) + (r + 0.5*sigma0^2)*T)/(sigma0*sqrt(T));
d2 = d1 - sigma0*sqrt(T);
N1 = normcdf(d1);
N2 = normcdf(d2);
European_call_BS = S0*N1 - K*exp(-r*T)*N2
%}

%function error=localvol