clear;

r1 = 0.06;
S01=17099.4;
A = importdata('Data_BNPP.txt','\t',1);
B1=A.data(:,:);

B1(:,2)=B1(:,2)/S01;
S01=1;
times1=unique(B1(:,1));
matur1=2;
a1=zeros(matur1,1);
b1=zeros(matur1,1);
c1=zeros(matur1,1);
d1=zeros(matur1,1);
sigmamax1=2;

M1 = 100; % number of asset paths
iterations1=2;

P1=B1;
for i=1:size(B1,1)
P1(i,3)=european_bs(S01,B1(i,2),r1,B1(i,3),B1(i,1)./252,'call');
end



%a1(1)=0.1963;b1(1)=17.9684;c1(1)=182.9531;
for iter=1:matur1
tic
ti1=times1(iter);
T1 = ti1/252;
D1=252;
L1 = T1*D1*2;

B1tmp=P1(P1(:,1)==times1(iter),2:3);
fun = @(var)Localvol(var(1),var(2),var(3),var(4),iter,S01,r1,T1,D1,M1,L1,B1tmp,a1,b1,c1,d1,sigmamax1,iterations1,times1);

A = [];b = [];Aeq = [];beq = [];nonlcon=[];
lb = [0.19,0.05,0,0.9];
ub = [0.25,0.6,0.5,1.2];
x0=[0.24-0.002*ti1, 0.2+0.3/ti1, 0.5/iter, 1.15];

options = optimoptions('patternsearch','Display','off','MaxIter',10000,'UseParallel',true);
%options = optimset('MaxFunEvals',100000,'MaxIter',100000,'Display','off');
vars=patternsearch(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
disp(vars);
error=Localvol(var(1),var(2),var(3),var(4),iter,S01,r1,T1,D1,M1,L1,B1tmp,a1,b1,c1,d1,sigmamax1,iterations1,times1);
a1(iter)=vars(1);
b1(iter)=vars(2);
c1(iter)=vars(3);
d1(iter)=vars(4);

time=toc;
if time>60
    time=strcat(num2str(floor(time/60)),"min");
else
    time=strcat(num2str(floor(time)),"sec");
end

format shortg;
c = clock;
disp(strcat("error=",num2str(error),strcat("   ",strcat(strcat(time,strcat("   ",strcat(num2str(c(4),'%02.f'),strcat(":",num2str(c(5),'%02.f')))))))))
fprintf('____________________________________________________\n\n')
end


sigma0=0.2;
iterations2=iterations1;
figure

for iter=1:matur1
ax(iter) = subplot(2,ceil(matur1/2),iter);
ti=times1(iter);
C=B1(B1(:,1)==ti,2:3);

T1 = times1(iter)/252;
D1=252;
L1 = T1*D1*2;

for i=1:size(C,1)
Euro(i)=Pricer_Vol(a1,b1,c1,d1,sigmamax1,S01,r1,T1,D1,M1,L1,C(i,1),iterations2,times1);
end

for i=1:size(C,1)
Euro_Const(i)=Pricer_Const(sigma0,S01,r1,T1,M1,L1,C(i,1),iterations2);
end


scatter(ax(iter),C(:,1),C(:,2),'.');
hold on;
scatter(ax(iter),C(:,1),Euro(:),'x');
hold on;
%scatter(ax(iter),C(:,1),Euro_Const(:));
title(ax(iter),times1(iter)) 
clear Euro Euro_Const
end

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


function Error_final=Localvol(ai,bi,ci,di,matur,S0,r,T,D,M,L,B,a,b,c,d,sigmamax,iterations,times)

dt = T/L;
Error=zeros(iterations,1);

a(matur)=ai;
b(matur)=bi;
c(matur)=ci;
d(matur)=di;

parfor iter=1:iterations
S = S0*ones(M,1); % asset paths
sigma=a(1)*ones(M,1);

for k = 2:L+1
S(:)=S(:)+S(:)*r*dt+sqrt(dt)*sigma(:).*S(:).*randn(M,1);

for i=1:size(times,1)
if dt*D*(k-1)<=times(i)
   %sigma=arrayfun(@(x) min(a(i)+b(i)*min(x,0)^2,sigmamax),(S(:,k)-c(i))/c(i));
   %sigma=arrayfun(@(x) min(a(i)+b(i)*(x-c(i))^2,sigmamax),(S(:)));
   sigma=arrayfun(@(x) min(a(i)+(-b(i)*(abs(x-d(i))-(x-d(i)))+c(i)*(abs(x-d(i))+(x-d(i)))).*(x-d(i)),sigmamax),(S(:)));
   %sigma=arrayfun(@(x) a(i)+b(i)*x^2,(S(:,k)-c(i))/c(i));
   break;
end
end
end
   
   Y=zeros(M,size(B,1));
   for j=1:size(B,1)
       for i=1:M
           Y(i,j) = max(S(i)-B(j,1),0);
       end
   end

   Error(iter)=sum((exp(-r*T)*mean(Y)-B(:,2)').^2)/(size(B,1));
   
end
Error_final=mean(Error);
end




function Euro_final=Pricer(a,b,c,d,sigmamax,S0,r,T,D,M,L,K,iterations,times)
dt = T/L;

for iter=1:iterations
S = S0*ones(M,1); % asset paths
sigma=a(1)*ones(M,1);

for k = 2:L+1
S(:)=S(:)+S(:)*r*dt+sqrt(dt)*sigma(:).*S(:).*randn(M,1);

for i=1:size(times,1)
if dt*D*(k-1)<=times(i)
   %sigma=arrayfun(@(x) min(a(i)+b(i)*min(x,0)^2,sigmamax),(S(:,k)-c(i))/c(i));
   sigma=arrayfun(@(x) min(a(i)+(-b(i)*(abs(x-d(i))-(x-d(i)))+c(i)*(abs(x-d(i))+(x-d(i)))).*(x-d(i)),sigmamax),(S(:)));
   %sigma=arrayfun(@(x) a(i)+b(i)*x^2,(S(:,k)-c(i))/c(i));
   break;
end
end

end

for i=1:M
Y(i) = max(S(i)-K,0);
end

Euro_Vol(iter)=exp(-r*dt*L)*mean(Y(:));

end

Euro_final=mean(Euro_Vol);
end


function Euro_final=Pricer_Vol(a,b,c,d,sigmamax,S0,r,T,D,M,L,K,iterations,times)
dt = T/L;

for iter=1:iterations
S = S0*ones(M,1); % asset paths
sigma=a(1)*ones(M,1);

for k = 2:L+1
S(:)=S(:)+S(:)*r*dt+sqrt(dt)*sigma(:).*S(:).*randn(M,1);

for i=1:size(times,1)
if dt*D*(k-1)<=times(i)
   %sigma=arrayfun(@(x) min(a(i)+b(i)*min(x,0)^2,sigmamax),(S(:,k)-c(i))/c(i));
   sigma=arrayfun(@(x) min(a(i)+(-b(i)*(abs(x-d(i))-(x-d(i)))+c(i)*(abs(x-d(i))+(x-d(i)))).*(x-d(i)),sigmamax),(S(:)));
   %sigma=arrayfun(@(x) a(i)+b(i)*x^2,(S(:,k)-c(i))/c(i));
   break;
end
end

end

for i=1:M
Y(i) = max(S(i)-K,0);
end

   euro=@(sigma)european_bs(S0,K,r,sigma,T,'call')-exp(-r*T)*mean(Y(:));
   Euro_Vol(iter)=fzero(euro,0.2);

end

Euro_final=mean(Euro_Vol);
end


function Euro_final_Const=Pricer_Const(sigma,S0,r,T,M,L,K,iterations)
dt = T/L;

for iter=1:iterations
S = S0*ones(M,L+1); % asset paths

for k = 2:L+1
S(:,k)=S(:,k-1)+S(:,k-1)*r*dt+sqrt(dt)*sigma.*S(:,k-1).*randn(M,1);

end

for i=1:M
Y(i) = max(S(i,L+1)-K,0);
end

Euro_Vol(iter)=exp(-r*dt*L)*mean(Y(:));
end

Euro_final_Const=mean(Euro_Vol);
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