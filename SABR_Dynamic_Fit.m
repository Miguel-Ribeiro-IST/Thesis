%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   FOR COMMENTS ON THE CODE, CHECK SIMILAR FILE SABR_Static_Fit.m   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
A = importdata('Data_BNPP.txt','\t',1);
B=A.data(:,:);


S0=17099.4;
r = 0.06;
matur=4;           %maturity until which we want to fit the data.
                   %If matur=5, all maturities until the fifth maturity in the file are chosen.
M=50000;
iterations=10;
beta=0.5;


B(:,2)=B(:,2)/S0;
S0=1;
B(:,1)=B(:,1)/252;
times=unique(B(:,1));
B=B(B(:,1)<=times(matur),:);

x0=[0.2,-0.1,0.1,1,0.1];

optimvars=Optimizer(beta,S0,B,r,x0);
alpha=optimvars(1);
rho0=optimvars(2);
nu0=optimvars(3);
a=optimvars(4);
b=optimvars(5);


Plotter(alpha,rho0,nu0,a,b,beta,S0,r,B,M,iterations,matur)



function optimvars=Optimizer(beta,S0,B,r,x0)
fun=@(var)SABRcal(var(1),var(2),var(3),var(4),var(5),beta,S0,B,r);
lb = [0,-1,0,0,0];
ub = [Inf,1,Inf,Inf,Inf];

options = optimoptions('patternsearch','Display','off');
optimvars=patternsearch(fun,x0,[],[],[],[],lb,ub,[],options);

fprintf(strcat("alpha=",strcat(num2str(optimvars(1)),strcat(",    rho0=",strcat(num2str(optimvars(2)),strcat(",    nu0=",strcat(num2str(optimvars(3)),strcat(",    a=",strcat(num2str(optimvars(4)),strcat(",    b=",strcat(num2str(optimvars(5)),strcat("\nerror=",strcat(num2str(SABRcal(optimvars(1),optimvars(2),optimvars(3),optimvars(4),optimvars(5),beta,S0,B,r)),"\n")))))))))))));
end


function Plotter(alpha,rho0,nu0,a,b,beta,S0,r,B,M,iterations,matur)
figure
times=unique(B(:,1));
for iter=1:matur
    ax(iter) = subplot(2,ceil(matur/2),iter);
    T=times(iter);
    
    C=B(B(:,1)==T,2:3);
    
    SABRVol=@(K)sigmaSABR(alpha,rho0,nu0,a,b,beta,K,S0*exp(r*T),T);
    
    for i=1:size(C,1)
        Volatility(i)= Pricer(alpha,rho0,nu0,a,b,beta,S0,r,T,M,T*252*2,C(i,1),iterations,"vol");
    end
    
    scatter(ax(iter),C(:,1),C(:,2),'.');
    hold on;
    scatter(ax(iter),C(:,1),Volatility(:),'x');
    hold on;
    fplot(ax(iter),SABRVol,[min(C(:,1)) max(C(:,1))])
    hold on;
    title(ax(iter),strcat(strcat(strcat(num2str(T*252)," days  ("),num2str(T*252/21))," months)"))
    clear Volatility
end
text1=strcat(strcat(strcat(num2str(times(1)*252)," days  ("),num2str(times(1)*252/21))," months)");
vars1=strcat(strcat(strcat("\beta=",num2str(beta)),strcat(",  paths=",num2str(M))),strcat(",  iterations=",num2str(iterations)));
text2=strcat(strcat(strcat(num2str(times(2)*252)," days  ("),num2str(times(2)*252/21))," months)");
vars2=strcat(strcat(strcat(strcat(strcat(strcat(strcat("\alpha=",num2str(alpha)),strcat(",  \rho0=",num2str(rho0))),strcat(",  \nu0=",num2str(nu0))),",  a="),num2str(a)),",  b="),num2str(b));
title(ax(1),{vars1,text1})
title(ax(2),{vars2,text2})
end


function Total_Error=SABRcal(alpha,rho,nu,a,b,beta,S0,B,r)
LS=0;
for i=1:size(B,1)
    LS=LS+((B(i,3)-sigmaSABR(alpha,rho,nu,a,b,beta,B(i,2),S0*exp(r*B(i,1)),B(i,1))))^2;
end
Total_Error=LS;
end



function sigma=sigmaSABR(alpha,rho0,nu0,a,b,beta,K,f,T)
w=alpha^(-1)*f^(1-beta);
n1=@(T)2*nu0*rho0/(T^2*(a+b)^2)*(exp(-(a+b)*T)-(1-(a+b)*T));
n22=@(T)3*nu0^2*rho0^2/(T^4*(a+b)^4)*(exp(-2*(a+b)*T)-8*exp(-(a+b)*T)+(7+2*(a+b)*T*(-3+(a+b)*T)));
v12=@(T)6*nu0^2/(2*b*T)^3*(((2*b*T)^2/2-2*b*T+1)-exp(-2*b*T));
v22=@(T)6*nu0^2/(2*b*T)^3*(2*(exp(-2*b*T)-1)+2*b*T*(exp(-2*b*T)+1));
%rhot=@(t)rho*exp(-a*t);
%nut=@(t)nu*exp(-b*t);

A1=@(T)(beta-1)/2+n1(T)*w/2;
A2=@(T)(1-beta)^2/12+(1-beta-n1(T)*w)/4+(4*v12(T)+3*(n22(T)-3*(n1(T))^2))*w^2/24;
B=@(T)1/w^2*((1-beta)^2/24+w*beta*n1(T)/4+(2*v22(T)-3*n22(T))*w^2/24);

sigma=1/w*(1+A1(T)*log(K./f)+A2(T)*(log(K./f)).^2+B(T)*T);
end


function Result_Avg=Pricer(alpha,rho0,nu0,a,b,beta,S0,r,T,M,L,K,iterations,PriceVol)
dt = T/L;
f=S0*exp(r*T);

parfor iter=1:iterations
    F = f*ones(M,1);
    alp=alpha*ones(M,1);
    
    for k = 1:L
        rho=rho0*exp(-a*dt*(k-1));
        nu=nu0*exp(-b*dt*(k-1));
        Z1=randn(M,1);
        Z2=rho*Z1+sqrt(1-rho^2)*randn(M,1);
        v=alp.*(F.^(beta-1));
        F(:)=F(:).*exp(v.*Z2.*sqrt(dt)-v.^2*dt/2);
        alp(:)=alp(:).*exp(nu*sqrt(dt)*Z1-nu^2*dt/2);
    end
    
    Y=zeros(M,1);
    for i=1:M
        Y(i) = max(F(i)-K,0);
    end
    
    if PriceVol=="price"
        Result(iter)=exp(-r*T)*mean(Y(:));
    else
        volatility=@(sigma)european_bs(S0,K,r,sigma,T,"call")-exp(-r*T)*mean(Y(:));
        Result(iter)=fzero(volatility,0.25);
    end
end

Result_Avg=mean(Result);
end


function price=european_bs(S0,K,r,sigma,T,putcall)
d1 = (log(S0/K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
N1 = normcdf(d1);
N2 = normcdf(d2);
if putcall=="call"
    price = S0*N1 - K*exp(-r*T)*N2;
elseif putcall=="put"
    price = S0*N1 - K*exp(-r*T)*N2 + K*exp(-r*T) - S0;
end
end