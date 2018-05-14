clear; %clear all variables. This prevents multiple simulations from interfering with one another
%Import data from file "Data_BNPP.txt" and assign it to matrix B. The file must be in the same folder
%1st column - maturities, 2nd column - strikes, 3rd column - implied volatilities
A = importdata('Data_BNPP.txt','\t',1);
B=A.data(:,:);

%%%%%%%%%%%%%%%%%%%%  INPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%
S0=17099.4;        %initial stock price
r = 0.06;          %risk-free rate. Forward prices in data file assumed r=0.06
matur=4;           %maturity at which we want to fit the data. If matur=5, the fifth maturity in the file is chosen.
M=50000;           %number of paths to be simulated
iterations=10;     %number of repetitions to be simulated (and then averaged)



%%%%%%%%%%%%%      ORIGINAL DATA MODIFICATIONS     %%%%%%%%%%%%
B(:,2)=B(:,2)/S0;     %normalize strike prices
S0=1;
B(:,1)=B(:,1)/252;    %convert maturities from days to years
times=unique(B(:,1));
B=B(B(:,2)>=0.7 ,:);     %only keep values of the maturity
%T=times(matur);       %selected maturity
%B=B(B(:,1)==T,:);     %only keep values of the maturity
%L = T*252*2;          %number of steps in simulations

P=B;
for i=1:size(B,1)
    P(i,3)=european_bs(S0,B(i,2),r,B(i,3),B(i,1),"call");
end


%%%%%%%%%%%%%%%%%%%%%    CALIBRATION      %%%%%%%%%%%%%%%%%%%%%%
x0=[5,0.05,0.08,-0.35,1.5];   %parameter starting values [kappa,nubar,nu0,rho,chi]

tic
optimvars=Optimizer(B,S0,r,x0);
toc

kappa=optimvars(1);
nubar=optimvars(2);
nu0=optimvars(3);
rho=optimvars(4);
chi=optimvars(5);



Plotter(kappa,nubar,nu0,rho,chi,S0,r,B,P,M,iterations,matur)

%HestonPrice(1.1,0.5,S0,r,kappa,nubar,nu0,rho,chi)



%lb = [0,0,0,-1,0];       %parameter lower bounds
%ub = [Inf,Inf,Inf,1,Inf];    %parameter upper bounds

%optimvars=Optimizer(beta,S0,B,r,x0);  %Obtain optimization variables

%Plot optimization results
%Plotter(alpha,rho,nu,beta,S0,r,T,L,M,iterations,B)


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

function result=HestonPrice(K,T,S0,r,kappa,nubar,nu0,rho,chi,PriceVol)
warning('off','all');
F=S0.*exp(r.*T);
fun1=@(u)real(exp(-1i.*u.*log(K))./(1i.*u.*F).*CharFuncHeston(u-1i,T,S0,r,kappa,nubar,nu0,rho,chi));
fun2=@(u)real(exp(-1i.*u.*log(K))./(1i.*u).*CharFuncHeston(u,T,S0,r,kappa,nubar,nu0,rho,chi));
P1=1/2+1/pi.*integral(fun1,0,100);
P2=1/2+1/pi.*integral(fun2,0,100);
call=S0.*P1-exp(-r.*T).*K.*P2;
if PriceVol=="price"
    result=call;
else
    volatility=@(sigma)european_bs(S0,K,r,sigma,T,"call")-call;
    result=fzero(volatility,0.25);  %Calculate the expected implied volatility
end
end

function w=CharFuncHeston(u,T,S0,r,kappa,nubar,nu0,rho,chi)
xi=kappa-chi.*rho.*1i.*u;
d=sqrt(xi.^2+chi.^2.*(u.^2+1i.*u));
A1=(u.^2+1i.*u).*sinh(d.*T/2);
%A2=d*cosh(d*T/2)+(kappa-chi*rho*1i*u)*sinh(d*t/2);
A2=exp(d.*T/2).*((d+xi)./2+(d-xi)./2.*exp(-d.*T));
A=A1./A2;
D=log(d)+kappa.*T/2-log(A2);
w=exp(1i.*u.*(log(S0)+r.*T)-T.*kappa.*nubar.*rho.*1i.*u./chi-nu0.*A+2*kappa.*nubar./chi.^2.*D);
end


%%%%%%%%%%%%%      DEFINE MINIMIZATION PROCEDURE      %%%%%%%%%%%%%%%
function optimvars=Optimizer(P,S0,r,x0)
warning('off','all');
fun=@(var)Hestoncal(var(1),var(2),var(3),var(4),var(5),S0,P,r); %function to be optimized.
%Hestoncal(kappa,nubar,nu0,rho,chi,S0,P,r)
lb = [0.0001,0.0001,0.0001,-1,0.0001];       %parameter lower bounds
ub = [100,5,5,1,5];    %parameter upper bounds


%options = optimoptions('patternsearch','Display','off'); %procedure options
%optimvars=patternsearch(fun,x0,[],[],[],[],lb,ub,@Feller_Condition,options); %define minimization procedure (patternsearch)
%optimvars=patternsearch(fun,x0,[],[],[],[],lb,ub,options); %define minimization procedure (patternsearch)

%options = optimoptions('simulannealbnd','Display','off'); %procedure options
%optimvars=simulannealbnd(fun,x0,lb,ub,options); %define minimization procedure (patternsearch)

rng default % For reproducibility
opts = optimoptions(@fmincon,'Display','off','Algorithm','sqp');
%problem = createOptimProblem('fmincon','objective',fun,'x0',x0,'lb',lb,'ub',ub,'options',opts,'nonlcon',@Feller_Condition);
problem = createOptimProblem('fmincon','objective',fun,'x0',x0,'lb',lb,'ub',ub,'options',opts);

ms = MultiStart('UseParallel',true,'StartPointsToRun','bounds-ineqs','Display','off');
[optimvars,f] = run(ms,problem,100);
%ms = GlobalSearch('StartPointsToRun','bounds-ineqs','Display','off');
%[optimvars,f] = run(ms,problem);

%print optimization output
fprintf(strcat("kappa=",strcat(num2str(optimvars(1)),strcat(",    nubar=",strcat(num2str(optimvars(2)),strcat(",    nu0=",strcat(num2str(optimvars(3)),strcat(strcat(strcat(",    rho=",num2str(optimvars(4))),(strcat(",    chi=",num2str(optimvars(5))))),strcat("\nerror=",strcat(num2str(f),"\n"))))))))));
end

function [c, ceq] = Feller_Condition(var)
c = var(5)^2-2*var(1)*var(2); %chi^2-2*kappa*nubar<=0
ceq = [];
end



%%%  CALCULATE LEAST SQUARES ERROR BETWEEN MODEL AND DATA  IMPLIED VOL  %%%
function Total_Error=Hestoncal(kappa,nubar,nu0,rho,chi,S0,P,r)
LS=0;
for i=1:size(P,1)
    %LS is the least squares error
    %Function sigmaSABR outputs the error between model and data implied
    %volatilities using the original Hagan formula
    %LS=LS+((P(i,3)-HestonPrice(P(i,2),P(i,1),S0,r,kappa,nubar,nu0,rho,chi,"price"))^2);
    LS=LS+(((P(i,3)-HestonPrice(P(i,2),P(i,1),S0,r,kappa,nubar,nu0,rho,chi,"vol"))./P(i,3))^2);
    %sigmaSABR(alpha,rho,nu,beta,K,f,T)
end
Total_Error=LS;
end




%%%%%%%%%%%%%%%%%%%%%    PLOT OPTIMIZATION RESULTS    %%%%%%%%%%%%%%%%%%%%
function Plotter(kappa,nubar,nu0,rho,chi,S0,r,B,P,M,iterations,matur)
figure
times=unique(B(:,1));
for iter=1:matur
    ax(iter) = subplot(2,ceil(matur/2),iter);
    T=times(iter);
    
    C=B(B(:,1)==T,2:3);
    %C=P(P(:,1)==T,2:3);
    
    %HestonP=@(K)HestonPrice(K,T,S0,r,kappa,nubar,nu0,rho,chi,"price");
    HestonVol=@(K)HestonPrice(K,T,S0,r,kappa,nubar,nu0,rho,chi,"vol");
    
    for i=1:size(C,1)
        %Price(i)= Pricer(kappa,nubar,nu0,rho,chi,C(i,1),S0,r,T,T*252*10,M,iterations,"price");
        Volatility(i)= Pricer(kappa,nubar,nu0,rho,chi,C(i,1),S0,r,T,T*252*10,M,iterations,"vol");
    end
    
    scatter(ax(iter),C(:,1),C(:,2),'.');
    hold on;
    scatter(ax(iter),C(:,1),Volatility(:),'x');
    %scatter(ax(iter),C(:,1),Price(:),'x');
    hold on;
    fplot(ax(iter),HestonVol,[0.5,1.5])
    %fplot(ax(iter),HestonP,[0.9,1.1])
    hold on;
    title(ax(iter),strcat(strcat(strcat(num2str(T*252)," days  ("),num2str(T*252/21))," months)"))
    clear Volatility
end
text1=strcat(strcat(strcat(num2str(times(1)*252)," days  ("),num2str(times(1)*252/21))," months)");
vars1=strcat(strcat("paths=",num2str(M)),strcat(",  iterations=",num2str(iterations)));
text2=strcat(strcat(strcat(num2str(times(2)*252)," days  ("),num2str(times(2)*252/21))," months)");
vars2=strcat("\kappa=",strcat(num2str(kappa),strcat(", \nu_{mean}=",strcat(num2str(nubar),strcat(", \nu_0=",strcat(num2str(nu0),strcat(strcat(strcat(", \rho=",num2str(rho)),(strcat(", \chi=",num2str(chi)))))))))));
title(ax(1),{vars1,text1})
title(ax(2),{vars2,text2})

end



%%% CALCULATE MONTE CARLO PRICE/IMPLIED VOLATILITY OF EUROPEAN OPTION UNDER SABR %%%
%Options are assumed Call
%If output should be a price, PriceVol="price"
%If output should be an implied volatility, PriceVol="vol"
function Result_Avg=Pricer(kappa,nubar,nu0,rho,chi,K,S0,r,T,L,M,iterations,PriceVol)
dt = T/L;      %time steps

%%NOTE:
%%%%We don't need to simulate an entire matrix of forwards for all time steps and for all paths.
%%%%We just need to simulate one vector of forwards for all paths and update it at each time step
parfor iter=1:iterations    %Perform the "for" cycle in parallel
    x = zeros(M,1);        %define initial vector of forwards
    nu=nu0*ones(M,1);    %define the initial vector of volatilities
    
    for k = 1:L
        Z1=randn(M,1);                                 %vector of random variables
        Z2=rho*Z1+sqrt(1-rho^2)*randn(M,1);            %vector of random variables with correlation "rho" with vector Z1
        x(:)=x(:)+(r-max(nu(:),0)./2)*dt+sqrt(dt.*max(nu(:),0)).*Z1;      %update forwards vector
        nu(:)=nu(:)+kappa.*(nubar-max(nu(:),0)).*dt+chi*sqrt(dt*max(nu(:),0)).*Z2+chi^2/4*dt*(Z2.^2-1);  %update volatilities vector
    end
    
    S=S0.*exp(x);
    
    Y=zeros(M,1);
    for i=1:M
        Y(i) = max(S(i)-K,0);  %Calculate the payoff of all paths (assuming calls)
    end
    
    if PriceVol=="price"
        Result(iter)=exp(-r*T)*mean(Y(:));  %Output the discounted expected payoff
    else
        volatility=@(sigma)european_bs(S0,K,r,sigma,T,"call")-exp(-r*T)*mean(Y(:));
        Result(iter)=fzero(volatility,0.25);  %Calculate the expected implied volatility
    end
    
end
Result_Avg=mean(Result);  %average all results
end


%%%%%%%%%%%%%%  CALCULATE BLACK-SCHOLES PRICE  %%%%%%%%%%%%%%%%%%%%
%If option is a call: putcall="call"
%If option is a put: putcall="put"
function price=european_bs(S0,K,r,sigma,T,putcall)
d1 = (log(S0./K) + (r + 0.5.*sigma.^2).*T)/(sigma.*sqrt(T));
d2 = d1 - sigma.*sqrt(T);
N1 = normcdf(d1);
N2 = normcdf(d2);
if putcall=="call"
    price = S0.*N1 - K.*exp(-r.*T).*N2;
elseif putcall=="put"
    price = S0.*N1 - K.*exp(-r*T).*N2 + K.*exp(-r.*T) - S0;
end
end