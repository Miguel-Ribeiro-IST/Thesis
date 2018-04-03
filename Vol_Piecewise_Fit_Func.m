clear;
tic
S01=167.65;
r1 = 0.06;


A = importdata('AAPL.txt','\t',1);
B1=A.data(:,:);

times1=unique(B1(:,1));

a1=zeros(size(times1,1),1);
b1=zeros(size(times1,1),1);
c1=zeros(size(times1,1),1);
sigmamax1=0.8;
matur1=5;
M1 = 3000; % number of asset paths
iterations1=150;
tic



%B1=B1(B1(:,1)==6,:);

%a1(1)=0.1963;b1(1)=17.9684;c1(1)=182.9531;
for iter=1:matur1

T1 = times1(iter)/252;
D1=252; %number of (trading) days in a year
L1 = T1*D1*2; % number of time intervals

fun = @(var)Localvol(var(1),var(2),var(3),iter,S01,r1,T1,D1,M1,L1,B1,sigmamax1,a1,b1,c1,iterations1);

A = [];b = [];Aeq = [];beq = [];nonlcon=[];
lb = [0.15,0.1,175];
ub = [0.3,25,190];
x0=[0.21 18/iter 183];

options = optimset('MaxFunEvals',100000,'MaxIter',100000,'Display','off','FinDiffRelStep',[0.001,0.01,0.01]);
vars=fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
a1(iter)=vars(1);
b1(iter)=vars(2);
c1(iter)=vars(3);
end

toc


ti=times1(matur1);
C=B1(B1(:,1)==ti,2:3);

for i=1:size(C,1)
Euro(i)=Pricer(a1,b1,c1,S01,r1,T1,D1,M1,L1,B1,sigmamax1,C(i,1),iterations1);
end

scatter(C(:,1),C(:,2));
hold on;
scatter(C(:,1),Euro(:));
beep
%{
for i=1:size(C,1)
euro=@(sigma)european_bs(S01,C(i,1),r1,sigma,ti/252,'call')-C(i,2);
C(i,2)=fzero(euro,0.5);
end
scatter(C(:,1),C(:,2));
hold on;
F = @(var,x)var(1)+var(2)*((x-var(3))/var(3)).^2;
data0 = [vars(1) vars(2) vars(3)];
opts=optimset('Display','off');
[var] = lsqcurvefit(F,data0,C(:,1),C(:,2),[0.1 0.1 165],[0.35 100 200],opts);
vol=@(x)var(1)+var(2)*((x-var(3))/var(3)).^2;
fplot(vol, [150 185]);
hold on;
vol2=@(x)vars(1)+vars(2)*((x-vars(3))/vars(3)).^2;
fplot(vol2, [150 185]);
%}


function Error_final=Localvol(at,bt,ct,matur,S0,r,T,D,M,L,B,sigmamax,a,b,c,iterations)

dt = T/L;
Error=zeros(iterations,1);
times=unique(B(:,1));

a(matur)=at;
b(matur)=bt;
c(matur)=ct;

for iter=1:iterations
S = S0*ones(M,L+1); % asset paths
sigma=a(1)*ones(M,1);

for k = 2:L+1
S(:,k)=S(:,k-1)+S(:,k-1)*r*dt+sqrt(dt)*sigma(:).*S(:,k-1).*randn(M,1);

for i=1:size(times,1)
if dt*D*(k-1)<times(i)
   sigma=arrayfun(@(x) min(a(i)+b(i)*max(x,0)^2,sigmamax),(S(:,k)-c(i))/c(i));
   break;
end
end

if dt*D*k==times(matur)%ismember(dt*D*k,times) 
   C=B(B(:,1)==dt*D*k,2:3);
   Y=zeros(M,size(C,1));
   for j=1:size(C,1)
       for i=1:M
           Y(i,j) = max(S(i,k)-C(j,1),0);
       end
   end
   Error(iter)=abs(sum(exp(-r*dt*k)*mean(Y)-C(:,2)'));
   clear Y;
end

end
end
Error_final=mean(Error);
end


function Euro_final=Pricer(a,b,c,S0,r,T,D,M,L,B,sigmamax,K,iterations)
dt = T/L;
times=unique(B(:,1));

for iter=1:iterations
S = S0*ones(M,L+1); % asset paths
sigma=a(1)*ones(M,1);

for k = 2:L+1
S(:,k)=S(:,k-1)+S(:,k-1)*r*dt+sqrt(dt)*sigma(:).*S(:,k-1).*randn(M,1);

for i=1:size(times,1)
if dt*D*(k-1)<times(i)
   sigma=arrayfun(@(x) min(a(i)+b(i)*x^2,sigmamax),(S(:,k)-c(i))/c(i));
   break;
end
end

end

for i=1:M
Y(i) = max(S(i,L+1)-K,0);
end

Euro_Vol(iter)=exp(-r*dt*L)*mean(Y(:));
end

Euro_final=mean(Euro_Vol);
end


function euro=european_bs(S0,K,r,sigma0,T,putcall)
d1 = (log(S0/K) + (r + 0.5*sigma0^2)*T)/(sigma0*sqrt(T));
d2 = d1 - sigma0*sqrt(T);
N1 = normcdf(d1);
N2 = normcdf(d2);
if putcall=='call'
    euro = S0*N1 - K*exp(-r*T)*N2;
elseif putcall=='put'
    euro = S0*N1 - K*exp(-r*T)*N2 + K*exp(-r*T) - S0;
end
end