clear;
S0=167.65;
r = 0.06;

T = 84/252;
D=252; %number of (trading) days in a year
L = T*D*2; % number of time intervals
M = 1000; % number of asset paths
sigmamax=1;

dt = T/L;
A = importdata('AAPL.txt','\t',1);
B=A.data(:,:);
Error=0;


times=unique(B(:,1));

%a=0.1897*ones(size(times,1),1);
a=zeros(size(times,1),1);
a(1)=0.1897;a(2)=0.1897;a(3)=0.1897;a(4)=0.1897;
%b=2.0986*ones(size(times,1),1);
b=zeros(size(times,1),1);
b(1)=20.0986;b(2)=20.0986;b(3)=20.0986;b(4)=20.0986;
c=179.6484*ones(size(times,1),1);


S = S0*ones(M,L+1); % asset paths

sigma=a(1)*ones(M,1);

%%{
for k = 2:L+1
S(:,k)=S(:,k-1)+S(:,k-1)*r*dt+sqrt(dt)*sigma(:).*S(:,k-1).*randn(M,1);

for i=1:size(times,1)
if dt*D*(k-1)<times(i)
   sigma=arrayfun(@(x) min(a(i)+b(i)*x^2,0.5),(S(:,k)-c(i))/c(i));
   break;
end
end

if ismember(dt*D*k,times) 
   C=B(B(:,1)==dt*D*k,2:3);
   for j=1:size(C,1)
       for i=1:M
           Y(i,j) = max(S(i,k)-C(j,1),0);
       end
   end
   Error=Error+abs(sum(exp(-r*dt*k)*mean(Y)-C(:,2)'));
   clear Y;
end

end
%}

plot(0:dt*D:T*D,S(1:100,:)')
Error
%{
Error=0;
% Find payoff Y at expiry.
for j=1:size(B,1)
for i=1:M
Y(i,j) = max(S(i,L+1)-B(j,2),0);
end
end
abs(sum(exp(-r*dt*L)*mean(Y)-B(:,3)'));
%}


%{
d1 = (log(S0/K) + (r + 0.5*sigma0^2)*T)/(sigma0*sqrt(T));
d2 = d1 - sigma0*sqrt(T);
N1 = normcdf(d1);
N2 = normcdf(d2);
European_call_BS = S0*N1 - K*exp(-r*T)*N2
%}

%function error=localvol