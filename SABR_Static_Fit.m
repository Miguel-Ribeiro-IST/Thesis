clear; %clear all variables. This prevents multiple simulations from interfering with one another
%Import data from file "Data_BNPP.txt" and assign it to matrix B1. The file must be in the same folder
%1st column - maturities, 2nd column - strikes, 3rd column - implied volatilities
A = importdata('Data_BNPP.txt','\t',1); 
B1=A.data(:,:);

%%%%%%%%%%%%%%%%%%%%  INPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%
r1 = 0.06;          %risk-free rate. Forward prices in data file assumed r=0.06
S01=17099.4;        %initial stock price
matur1=2;           %maturity at which we want to fit the data. If matur1=5, the fifth maturity in the file is chosen.
M1=50000;           %number of paths to be simulated
iterations=10;      %number of repetitions to be simulated (and then averaged)
beta=1;             %parameter of the SABR model


%%%%%%%%%%%%%      ORIGINAL DATA MODIFICATIONS     %%%%%%%%%%%%
B1(:,2)=B1(:,2)/S01;    %normalize strike prices
S01=1;
B1(:,1)=B1(:,1)/252;    %convert maturities from days to years
times1=unique(B1(:,1));
T1=times1(matur1);      %selected maturity
B1=B1(B1(:,1)==T1,:);   %only keep values of the maturity
L1 = T1*252*2;          %number of steps in simulations


%%%%%%%%%%%%%%%%%%%%%    CALIBRATION      %%%%%%%%%%%%%%%%%%%%%%
x0=[0.5,-0.5,0.5];   %parameter starting values [alpha, rho, nu]
optimvars=Optimizer(beta,S01,B1,r1,x0);  %Obtain optimization variables
alpha=optimvars(1);
rho=optimvars(2);
nu=optimvars(3);

%Plot optimization results
Plotter(alpha,rho,nu,beta,S01,r1,T1,L1,M1,iterations,B1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%       FUNCTIONS       %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   PARAMETERS   %%%%
%alpha, beta, rho, nu = SABR parameters
%S0=starting stock price
%f=starting forward price
%r=risk-free rate
%T=maturity (in years)
%K=strike price
%L=number of simulation time steps
%M=number of simulated paths
%iterations=number of repetitions of the simulations (to reduce error)
%B=input data file
%x0=optimization starting parameters

%%%%%%%%%%%%%      DEFINE MINIMIZATION PROCEDURE      %%%%%%%%%%%%%%%
function optimvars=Optimizer(beta,S0,B,r,x0)
fun=@(var)SABRcal(var(1),var(2),var(3),beta,S0,B,r); %function to be optimized.
%alpha=var(1), rho=var(2), nu=var(3);
%SABRcal(alpha,rho,nu,beta,S0,B,r)
lb = [0,-1,0];       %parameter lower bounds
ub = [Inf,1,Inf];    %parameter upper bounds

options = optimoptions('patternsearch','Display','off'); %procedure options
optimvars=patternsearch(fun,x0,[],[],[],[],lb,ub,[],options); %define minimization procedure (patternsearch)

%print optimization output
fprintf(strcat("alpha=",strcat(num2str(optimvars(1)),strcat(",    rho=",strcat(num2str(optimvars(2)),strcat(",    nu=",strcat(num2str(optimvars(3)),strcat("\nerror=",strcat(num2str(SABRcal(optimvars(1),optimvars(2),optimvars(3),beta,S0,B,r)),"\n")))))))));
end


%%%%%%%%%%%%%%%%%%%%%    PLOT OPTIMIZATION RESULTS    %%%%%%%%%%%%%%%%%%%%
function Plotter(alpha,rho,nu,beta,S0,r,T,L,M,iterations,B)
for i=1:size(B,1)
    P(i)=european_bs(S0,B(i,2),r,B(i,3),B(i,1),"call");  %obtain the BS price from each of the data's implied volatilities
    %european_bs(S0,K,r,sigma0,T,putcall)
    
    %Obtain Monte Carlo prices and implied volatilities
    Price(i)=Pricer(alpha,rho,nu,beta,B(i,2),S0*exp(r*T),r,T,L,M,iterations,"price");
    Volatility(i)=Pricer(alpha,rho,nu,beta,B(i,2),S0*exp(r*T),r,T,L,M,iterations,"vol");
    %Pricer(alpha,rho,nu,beta,K,f,r,T,L,M,iterations,PriceVol)
end

%Plot implied volatilities from data, from Monte-Carlo and the SABR implied volatility function
SABRVol=@(K)sigmaSABR(alpha,rho,nu,beta,K,S0*exp(r*T),T);   %SABR implied volatility as a function of K
ax1 = subplot(1,2,1);
scatter(ax1,B(:,2),B(:,3),'.');
hold on;
scatter(ax1,B(:,2),Volatility(:),'x');
hold on;
fplot(ax1,SABRVol,[min(B(:,2)) max(B(:,2))],'k')

%Plot prices from data and from Monte-Carlo
ax2 = subplot(1,2,2);
scatter(ax2,B(:,2),P(:),'.');
hold on;
scatter(ax2,B(:,2),Price(:),'x');

%Show fit results in plot titles
text1=strcat(strcat(strcat(strcat("\beta=",num2str(beta)),strcat(", paths=",num2str(M))),strcat(", iterations=",num2str(iterations))),strcat(", maturity=",num2str(T*252)));
text2=strcat(strcat(strcat("\alpha=",num2str(alpha)),strcat(", \rho=",num2str(rho))),strcat(", \nu=",num2str(nu)));
title(ax1,{text1,'Volatilities'})
title(ax2,{text2,'Prices'})
end


%%%  CALCULATE LEAST SQUARES ERROR BETWEEN MODEL AND DATA  IMPLIED VOL  %%%
function Total_Error=SABRcal(alpha,rho,nu,beta,S0,B,r)
LS=0;
for i=1:size(B,1)
    %LS is the least squares error
    %Function sigmaSABR outputs the error between model and data implied
    %volatilities using the original Hagan formula
    LS=LS+((B(i,3)-sigmaSABR(alpha,rho,nu,beta,B(i,2),S0*exp(r*B(i,1)),B(i,1)))^2);
    %sigmaSABR(alpha,rho,nu,beta,K,f,T)
end
Total_Error=LS;
end


%%%%%%%%%%% CALCULATE SABR MODEL IMPLIED VOLATILITY FROM HAGAN  %%%%%%%%%%
%%Find function in page 9 of file github.com/Miguel-Ribeiro-IST/Thesis/blob/master/References/Hagan_SABR.pdf
function sigma=sigmaSABR(alpha,rho,nu,beta,K,f,T)
z=nu./alpha.*(f.*K).^((1-beta)./2).*log(f./K);
x=log((sqrt(1-2.*rho.*z+z.^2)+z-rho)./(1-rho));
if f==K
    sigma=alpha./((f).^(1-beta)).*(1+T.*((1-beta)^2/24.*alpha^2/((f).^(2-2*beta))+1/4*rho.*beta.*alpha.*nu./((f).^(1-beta))+(2-3*rho^2)/24.*nu^2));
elseif beta==0
    sigma=alpha.*log(f./K)./(f-K).*(z./x).*(1+T.*(alpha.^2./(24.*f.*K)+(2-3*rho.^2)./24.*nu.^2));
elseif beta==1
    sigma=alpha.*(z./x).*(1+T.*(1/4*rho.*alpha.*nu+1/24*(2-3*rho.^2).*nu.^2));
else
    sigma=alpha./((f.*K).^((1-beta)./2).*(1+(1-beta).^2/24.*(log(f./K)).^2+(1-beta).^4./1920.*(log(f./K)).^4)).*(z./x).*(1+((1-beta).^2/24.*alpha.^2./(f.*K).^(1-beta)+1/4.*rho.*nu.*beta.*alpha./(f.*K).^((1-beta)./2)+(2-3*rho.^2)./24.*nu.^2).*T);
end
end


%%% CALCULATE MONTE CARLO PRICE/IMPLIED VOLATILITY OF EUROPEAN OPTION UNDER SABR %%%
%Options are assumed Call
%If output should be a price, PriceVol="price"
%If output should be an implied volatility, PriceVol="vol"
function Result_Avg=Pricer(alpha,rho,nu,beta,K,f,r,T,L,M,iterations,PriceVol)
dt = T/L;      %time steps

%%NOTE:
%%%%We don't need to simulate an entire matrix of forwards for all time steps and for all paths.
%%%%We just need to simulate a vector of forwards for all paths and update it at each time step
    
parfor iter=1:iterations    %Perform the "for" cycle in parallel
    F = f*ones(M,1);        %define initial vector of forwards
    alp=alpha*ones(M,1);    %define the initial vector of volatilities
    
    for k = 1:L
        Z1=randn(M,1);         %vector of random variables
        Z2=rho*Z1+sqrt(1-rho^2)*randn(M,1);     %vector of random variables with correlation "rho" with vector Z1
        v=alp.*(F.^(beta-1));                          %intermediate variable
        F(:)=F(:).*exp(v.*Z2*sqrt(dt)-v.^2*dt/2);      %update forwards vector
        alp(:)=alp(:).*exp(nu*sqrt(dt)*Z1-nu^2*dt/2);  %update volatilities vector
    end
    
    Y=zeros(M,1);
    for i=1:M
        Y(i) = max(F(i)-K,0);  %Calculate the payoff of all paths (assuming calls)
    end
    
    
    if PriceVol=="price"
        Result(iter)=exp(-r*T)*mean(Y(:));  %Output the discounted expected payoff
    else
        volatility=@(sigma)european_bs(f*exp(-r*T),K,r,sigma,T,"call")-exp(-r*T)*mean(Y(:));
        Result(iter)=fzero(volatility,0.25);  %Calculate the expectedimplied volatility
    end
    
end
Result_Avg=mean(Result);  %average all results
end


%%%%%%%%%%%%%%  CALCULATE BLACK-SCHOLES PRICE  %%%%%%%%%%%%%%%%%%%%
%If option is a call: putcall="call"
%If option is a put: putcall="put"
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