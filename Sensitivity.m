clear;
tic
S0 = 36;
K = 40;
T = 1;
r = 0.06;
sigma = 0.2;
L = 50; % number of time intervals
M = 50000; % number of asset paths
%for q=1:30
%val(q)=Price(S0,K,T,r,sigma,L,M);
%end
%[mean(val),std(val)]

N=1000;
Un=3;
rm=0.06;
rs=0.005;
sigmam=0.2;
sigmas=0.02;
S0m=36;
S0s=sigmam;
%rs=0;sigmas=0;S0s=0;

A=horzcat(K*ones(N,1),T*ones(N,1),L*ones(N,1),M*ones(N,1),normrnd(S0m,S0s,N,1),normrnd(rm,rs,N,1),normrnd(sigmam,sigmas,N,1));
fA=Price(A)';
B=horzcat(K*ones(N,1),T*ones(N,1),L*ones(N,1),M*ones(N,1),normrnd(S0m,S0s,N,1),normrnd(rm,rs,N,1),normrnd(sigmam,sigmas,N,1));
fB=Price(B)';


for j=1:Un
   C(:,:,j)=A;
   C(:,4+j,j)=B(:,4+j);
   fC(:,j)=Price(C(:,:,j))';
end

for i=1:Un
   Si(i)=1/N*sum(fB.*(fC(:,i)-fA));
   STi(i)=1/(2*N)*sum((fA-fC(:,i)).^2);
end

Si
STi
beep

toc