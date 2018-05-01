clear;
A = importdata('Data_BNPP.txt','\t',1);
B1=A.data(:,:);
S01=17099.4;
r1 = 0.06;


matur1=2;
M1 = 10000;
iterations1=50;
beta=1;
PriceVol="vol";
iterations2=iterations1;



B1(:,2)=B1(:,2)/S01;
S01=1;
B1(:,1)=B1(:,1)/252;
times1=unique(B1(:,1));
T1=times1(matur1);
B1=B1(B1(:,1)==T1,:);
L1 = T1*252*2;


P1=B1;
for i=1:size(B1,1)
    P1(i,3)=european_bs(S01,B1(i,2),r1,B1(i,3),B1(i,1),"call");
end


tic
if PriceVol=="price"
    B1tmp=P1(P1(:,1)==T1,2:3);
else
    B1tmp=B1(B1(:,1)==T1,2:3);
end

fun = @(var)SABRvol(var(1),var(2),var(3),beta,S01,r1,T1,M1,L1,B1tmp,iterations1,PriceVol);
lb = [0,-1,0];
ub = [Inf,1,Inf];
x0=[0.24603     -0.44043        1.476];
A = [];b = [];Aeq = [];beq = [];nonlcon=[];
options = optimoptions('patternsearch','Display','off','MaxIter',10000,'UseParallel',true);
vars=patternsearch(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
disp(vars);

alpha=vars(1);
rho=vars(2);
nu=vars(3);



Timer(SABRvol(var(1),var(2),var(3),beta,S01,r1,T1,M1,L1,B1tmp,iterations1,PriceVol),toc);



    C1=P1(:,2:3);
    C2=B1(:,2:3);

    SABRVol=@(K)sigmaSABR(alpha,rho,nu,beta,K,S01*exp(r1*T1),T1);

    
for i=1:size(C1,1)
    Euro(i)=Pricer(alpha,rho,nu,beta,S01,r1,T1,M1,L1,C1(i,1),iterations2,"price");
    Eurovol(i)=Pricer(alpha,rho,nu,beta,S01,r1,T1,M1,L1,C2(i,1),iterations2,"vol");
end

ax1 = subplot(1,2,1);
scatter(ax1,C1(:,1),C1(:,2),'.');
hold on;
scatter(ax1,C1(:,1),Euro(:),'x');

ax2 = subplot(1,2,2);
scatter(ax2,C2(:,1),C2(:,2),'.');
hold on;
scatter(ax2,C2(:,1),Eurovol(:),'x');
hold on;
fplot(ax2,SABRVol,[min(C1(:,1)) max(C1(:,1))])

text1=strcat(strcat(strcat(strcat("\beta=",num2str(beta)),strcat(", paths=",num2str(M1))),strcat(", iterations=",num2str(iterations1))),strcat(", maturity=",num2str(T1*252)));
text2=strcat(strcat(strcat("\alpha=",num2str(alpha)),strcat(", \rho=",num2str(rho))),strcat(", \nu=",num2str(nu)));
title(ax1,text1)
title(ax2,text2)



function Error_final=SABRvol(alpha,rho,nu,beta,S0,r,T,M,L,B,iterations,PriceVol)

dt = T/L;
Error=zeros(iterations,1);
f=S0*exp(r*T);

parfor iter=1:iterations
    F = f*ones(M,1);
    alp=alpha*ones(M,1);
    
    for k = 1:L
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
        Error(iter)=sum(((exp(-r*T)*mean(Y)-B(:,2)').^2)./((S0-B(:,2)').^2+0.25));
    else
        Z=exp(-r*T)*mean(Y);
        for i=1:size(B,1)
            euro=@(sigma)european_bs(S0,B(i,1),r,sigma,T,"call")-Z(i);
            Z(i)=fzero(euro,0.25);
        end
        Error(iter)=sum(((Z-B(:,2)').^2)./((S0-B(:,2)').^2+0.25));
    end
    
end
Error_final=mean(Error);
end



function Euro_final=Pricer(alpha,rho,nu,beta,S0,r,T,M,L,K,iterations,PriceVol)
dt = T/L;
f=S0*exp(r*T);
for iter=1:iterations
    F = f*ones(M,1);
    alp=alpha*ones(M,1);
    
    for k = 1:L
        Z1=randn(M,1);
        alp(:)=alp(:).*exp(nu*sqrt(dt)*Z1-nu^2*dt/2);
        v=alp.*(F.^(beta-1));
        F(:)=F(:).*exp(v.*(rho*Z1+sqrt(1-rho^2)*randn(M,1))*sqrt(dt)-v.^2*dt/2);
        
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


function sigma=sigmaSABR(alpha,rho,nu,beta,K,f,T)
z=nu./alpha.*(f.*K).^((1-beta)./2).*log(f./K);
x=log((sqrt(1-2.*rho.*z+z.^2)+z-rho)./(1-rho));
sigma=alpha./((f.*K).^((1-beta)./2).*(1+(1-beta).^2/24.*(log(f./K)).^2+(1-beta).^4./1920.*(log(f./K)).^4)).*(z./x).*(1+((1-beta).^2/24.*alpha.^2./(f.*K).^(1-beta)+1/4.*rho.*nu.*beta.*alpha./(f.*K).^((1-beta)./2)+(2-3*rho.^2)./24.*nu.^2).*T);
end