clear;
A = importdata('Data_BNPP.txt','\t',1);
B1=A.data(:,:);
S01=17099.4;
r1 = 0.06;


matur1=4;
M1 = 100;
iterations1=3;
beta=0.5;
PriceVol="vol";
iterations2=iterations1;
D1=2*252;



B1(:,2)=B1(:,2)/S01;
S01=1;
B1(:,1)=B1(:,1)/252;
times1=unique(B1(:,1));
B1=B1(B1(:,1)<=times1(matur1),:);
times1=unique(B1(:,1));

P1=B1;
for i=1:size(B1,1)
    P1(i,3)=european_bs(S01,B1(i,2),r1,B1(i,3),B1(i,1),"call");
end


tic
if PriceVol=="price"
    B1tmp=P1;
else
    B1tmp=B1;
end


tic
fun = @(var)DynSABRvol(var(1),var(2),var(3),var(4),var(5),beta,S01,r1,times1,M1,D1,B1tmp,iterations1,PriceVol);
lb = [0,-1,0,0,0];
ub = [Inf,1,Inf,Inf,Inf];
x0=[0.3,-0.5,0.1,1,1];
A = [];b = [];Aeq = [];beq = [];nonlcon=[];
options = optimoptions('patternsearch','Display','off','MaxIter',10000,'UseParallel',true);
vars=patternsearch(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
disp(vars);

alpha=vars(1);
rho=vars(2);
nu=vars(3);
a=vars(4);
b=vars(5);
Timer(DynSABRvol(alpha,rho,nu,a,b,beta,S01,r1,times1,M1,D1,B1tmp,iterations1,PriceVol),toc);
beep

figure
for iter=1:matur1
    ax(iter) = subplot(2,ceil(matur1/2),iter);
    ti=times1(iter);
    
    C=B1(B1(:,1)==ti,2:3);
    

    for i=1:size(C,1)
    Eurovol(i)=Pricer(alpha,rho,nu,a,b,beta,S01,r1,ti,M1,D1,C(i,1),iterations2,"vol");
end
    
    scatter(ax(iter),C(:,1),C(:,2),'.');
    hold on;
    scatter(ax(iter),C(:,1),Eurovol(:),'x');
    hold on;
    title(ax(iter),times1(iter)*252)
end


function Error_final=DynSABRvol(alpha,rho,nu,a,b,beta,S0,r,times,M,D,Btmp,iterations,PriceVol)
LS=0;
for i=1:size(times,1)
    C=Btmp(Btmp(:,1)==times(i),2:3);
    LS=LS+SABRvol(alpha,rho,nu,a,b,beta,S0,r,times(i),M,D,C,iterations,PriceVol);
end
Error_final=LS;
end



function Error_final=SABRvol(alpha,rho0,nu0,a,b,beta,S0,r,T,M,L,B,iterations,PriceVol)
dt = T/L;
Error=zeros(iterations,1);
f=S0*exp(r*T);

parfor iter=1:iterations
    F = f*ones(M,1);
    alp=alpha*ones(M,1);
    
    for k = 1:L
        rho=rho0*exp(-a*dt*k);
        nu=nu0*exp(-b*dt*k);
        Z1=randn(M,1);
        alp(:)=alp(:).*exp(nu*sqrt(dt)*Z1-nu^2*dt/2);
        v=alp.*(F.^(beta-1));
        F(:)=F(:).*exp(v.*(rho*Z1+sqrt(1-rho^2)*randn(M,1))*sqrt(dt)-v.^2*dt/2);
        
    end
    
    Y=zeros(M,size(B,1));
    for j=1:size(B,1)
        for i=1:M
            Y(i,j) = max(F(i)-B(j,1),0);
        end
    end
    
    if PriceVol=="price"
        Error(iter)=sum((exp(-r*T)*mean(Y)-B(:,2)').^2);
    else
        Z=exp(-r*T)*mean(Y);
        for i=1:size(B,1)
            euro=@(sigma)european_bs(S0,B(i,1),r,sigma,T,"call")-Z(i);
            Z(i)=fzero(euro,0.25);
        end
        Error(iter)=sum((Z-B(:,2)').^2);
    end
    
end
Error_final=mean(Error);
end



function Euro_final=Pricer(alpha,rho0,nu0,a,b,beta,S0,r,T,M,L,K,iterations,PriceVol)
dt = T/L;
f=S0*exp(r*T);

for iter=1:iterations
    F = f*ones(M,1);
    alp=alpha*ones(M,1);
    
    for k = 1:L
        rho=rho0*exp(-a*dt*k);
        nu=nu0*exp(-b*dt*k);
        Z1=randn(M,1);
        v=alp.*(F.^(beta-1));
        F(:)=F(:).*exp(v.*(rho*Z1+sqrt(1-rho^2)*randn(M,1))*sqrt(dt)-v.^2*dt/2);
        alp(:)=alp(:).*exp(nu*sqrt(dt)*Z1-nu^2*dt/2);
    end
    
    Y=zeros(M,1);
    for i=1:M
        Y(i) = max(F(i)-K,0);
    end
    
    if PriceVol=="price"
        Euro(iter)=exp(-r*T)*mean(Y(:));
    else
        euro=@(sigma)european_bs(S0,K,r,sigma,T,"call")-exp(-r*T)*mean(Y(:));
        Euro(iter)=fzero(euro,0.25);
    end
end

Euro_final=mean(Euro);
end



function euro=european_bs(S0,K,r,sigma0,T,putcall)
d1 = (log(S0/K) + (r + 0.5*sigma0^2)*T)/(sigma0*sqrt(T));
d2 = d1 - sigma0*sqrt(T);
N1 = normcdf(d1);
N2 = normcdf(d2);
if putcall=="call"
    euro = S0*N1 - K*exp(-r*T)*N2;
elseif putcall=="put"
    euro = S0*N1 - K*exp(-r*T)*N2 + K*exp(-r*T) - S0;
end
end

function Timer(error,time)
if time>3600
    time=strcat(strcat(strcat(num2str(floor(time/3600)),"hrs,"),num2str(floor((time-floor(time/3600)*3600)/60))),"min");
elseif time>60
    time=strcat(strcat(strcat(num2str(floor(time/60)),"min,"),num2str(floor(time-floor(time/60)*60))),"sec");
else
    time=strcat(num2str(floor(time)),"sec");
end
format shortg;
c = clock;
disp(strcat("error=",num2str(error),strcat("   ",strcat(strcat(time,strcat("   ",strcat(num2str(c(4),'%02.f'),strcat(":",num2str(c(5),'%02.f')))))))))
fprintf('____________________________________________________\n\n')
end