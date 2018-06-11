clear;


A = importdata('Data_BNPP.txt','\t',1);
B1=A.data(:,:);

r1 = 0.06;%%if edited, the data must be changed!
S01=17099.4;
matur1=2;
M1=50000;
iterations=50;
Mplot=100;
beta=0;

B1(:,2)=B1(:,2)/S01;
S01=1;
B1(:,1)=B1(:,1)/252;
times1=unique(B1(:,1));
B1=B1(B1(:,1)==times1(matur1),:);


P1=B1;
for i=1:size(B1,1)
    P1(i,3)=european_bs(S01,B1(i,2),r1,B1(i,3),B1(i,1),"call");
    %european_bs(S0,K,r,sigma0,T,putcall)
end


fun=@(var)SABRcal(var(1),var(2),var(3),beta,S01,B1,r1);
%SABRcal(alpha,rho,nu,beta,S0,B,r)
lb = [0,-1,0];
ub = [Inf,1,Inf];
x0=[1,0.5,1];
% fun=@(var)SABRcal(var(1),var(2),var(3),var(4),S01,B1,r1);
% lb = [0,-1,0,0];
% ub = [Inf,1,Inf,1];
% x0=[0.5,0.5,0.5,0.5];

A = [];b = [];Aeq = [];beq = [];nonlcon=[];
options = optimoptions('patternsearch','Display','off','MaxIter',10000,'UseParallel',true);
vars=patternsearch(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
disp(vars);
error=SABRcal(vars(1),vars(2),vars(3),beta,S01,B1,r1);
% error=SABRcal(vars(1),vars(2),vars(3),var(4),S01,B1,r1);
% beta=var(4);
disp(strcat("      error=",num2str(error)));

alpha=vars(1);
rho=vars(2);
nu=vars(3);


ti=times1(matur1);
C=B1(B1(:,1)==ti,2:3);
T1 = ti;
L1 = T1*252*2;


SABRVol=@(K)sigmaSABR(alpha,rho,nu,beta,K,S01*exp(r1*ti),ti);
C=P1(P1(:,1)==ti,2:3);
C2=B1(B1(:,1)==ti,2:3);

for i=1:size(C,1)
    Euro(i)=Pricer(alpha,rho,nu,beta,C(i,1),S01*exp(r1*ti),r1,T1,L1,M1,iterations,"price");
    EuroVol(i)=Pricer(alpha,rho,nu,beta,C(i,1),S01*exp(r1*ti),r1,T1,L1,M1,iterations,"vol");
end

ax1 = subplot(1,2,1);
scatter(ax1,C2(:,1),C2(:,2),'.');
hold on;
%scatter(C(:,1),SABRVol(:),'x');
fplot(ax1,SABRVol,[min(C2(:,1)) max(C2(:,1))])
hold on;
scatter(ax1,C2(:,1),EuroVol(:),'x');
%scatter(ax(iter),C(:,1),Euro_Const(:));
%}


ax2 = subplot(1,2,2);
scatter(ax2,C(:,1),C(:,2),'.');
hold on;
scatter(ax2,C(:,1),Euro(:),'x');
hold on;

%ax3 = subplot(2,2,2);
%Plotter(alpha,rho,nu,beta,S01*exp(r1*ti),r1,T1,L1,Mplot,ax3)

text1=strcat(strcat(strcat(strcat("\beta=",num2str(beta)),strcat(", paths=",num2str(M1))),strcat(", iterations=",num2str(iterations))),strcat(", maturity=",num2str(T1*252)));
text2=strcat(strcat(strcat("\alpha=",num2str(alpha)),strcat(", \rho=",num2str(rho))),strcat(", \nu=",num2str(nu)));
title(ax1,text1)
title(ax2,text2)


function Error=SABRcal(alpha,rho,nu,beta,S0,B,r)
LS=0;
for i=1:size(B,1)
    LS=LS+((B(i,3)-sigmaSABR(alpha,rho,nu,beta,B(i,2),S0*exp(r*B(i,1)),B(i,1)))/(B(i,3)))^2;
end
Error=LS;
end



function sigma=sigmaSABR(alpha,rho,nu,beta,K,f,T)
z=nu./alpha.*(f.*K).^((1-beta)./2).*log(f./K);
%z=nu.*(f.^(1-beta)-K.^(1-beta))./(alpha.*(1-beta));
x=log((sqrt(1-2.*rho.*z+z.^2)+z-rho)./(1-rho));
sigma=alpha./((f.*K).^((1-beta)./2).*(1+(1-beta).^2/24.*(log(f./K)).^2+(1-beta).^4./1920.*(log(f./K)).^4)).*(z./x).*(1+((1-beta).^2/24.*alpha.^2./(f.*K).^(1-beta)+1/4.*rho.*nu.*beta.*alpha./(f.*K).^((1-beta)./2)+(2-3*rho.^2)./24.*nu.^2).*T);
%sigma=1./(1+(1-beta).^2.*(log(f./K)).^2/24+(1-beta).^4.*(log(f./K)).^4/1920).*(nu.*log(f./K)./x).*(1+T*((1-beta).^2/24*alpha^2./(K.*f).^(1-beta)+1/4*rho.*beta.*nu.*alpha./(K.*f).^((1-beta)./2)+(2-3*rho.^2)/24.*nu.^2));
end


function Euro_final=Pricer(alpha,rho,nu,beta,K,f,r,T,L,M,iterations,PriceVol)
dt = T/L;
parfor iter=1:iterations
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
        euro=@(sigma)european_bs(f*exp(-r*T),K,r,sigma,T,"call")-exp(-r*T)*mean(Y(:));
        Euro(iter)=fzero(euro,0.25);
    end
    
end
Euro_final=mean(Euro);
end


function Plotter(alpha,rho,nu,beta,f,r,T,L,M,ax)
dt = T/L;
F = f*ones(M,L);
alp=alpha*ones(M,1);

for k = 1:L
    Z1=randn(M,1);
    alp(:)=alp(:).*exp(nu*sqrt(dt)*Z1-nu^2*dt/2);
    v=alp.*(F(:,k).^(beta-1));
    F(:,k+1)=F(:,k).*exp(v.*(rho*Z1+sqrt(1-rho)*randn(M,1))*sqrt(dt)-v.^2*dt/2);
end

for k = 1:(L+1)
    F(:,k)=F(:,k).*exp(-r*(T-dt*(k-1)));
end

plot(ax,0:dt*252:T*252,F(:,:)')
xlim([0 T*252])
ylim([0 2])

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