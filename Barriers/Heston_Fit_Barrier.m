clearvars -except kappa nubar nu0 rho chi; %clear all variables. This prevents multiple simulations from interfering with one another
%Import data from file "Data_BNPP_2.txt" and assign it to matrix B. The file must be in the same folder
%1st column - maturities, 2nd column - strikes, 3rd column - implied volatilities
kappa=53.4355;nubar=0.0653;nu0=0.1046;rho=-0.4086;chi=6.2554;
A = importdata('Data_BNPP.txt','\t',1);
B=A.data(:,:);

%%%%%%%%%%%%%%%%%%%%  INPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%
S0=17099.4;               %initial stock price
r = 0;                    %risk-free rate. Forward prices in data file assumed r=0.0
matur=2;                  %maturity until which we want to fit the data. If matur=5, all maturities until the fifth are chosen.
OptAlg="CMA";      %"CMA" or "MultiStart" optimization algorithms


%%%%%%%%%%%%%%%%%%%   MONTE CARLO SIMULATION %%%%%%%%%%%%%%
%After calibrating all the model's parameters, we may want to simulate the implied volatilities using Monte Carlo
SimPoints=true;   %true or false - define if Monte Carlo simulation should be executed
M=100000;           %number of paths to be simulated
repetitions=10;
barr=[1.05,1.2,1.4];
%L=T*252*2


%%%%%%%%%%%%%      ORIGINAL DATA MODIFICATIONS     %%%%%%%%%%%%
B(:,2)=B(:,2)/S0;                %normalize strike prices
S0=1;
B(:,1)=B(:,1)/252;               %convert maturities from days to years
times=unique(B(:,1));
T=times(matur);
B=B(B(:,1)==T,:);     %only keep values of the maturity


%%%%%%%%%%%%%%%%%%%%%    CALIBRATION      %%%%%%%%%%%%%%%%%%%%%%
%    [kappa, nubar, nu0, rho, chi]
x0=  [40.5,  0.065, 0.1, -0.4, 5];      %parameter starting values
lb = [0,     0, 0, -1,  0];          %parameter lower bounds
ub = [60,    2, 2, 1,   10];          %parameter upper bounds
%{
optimvars=Optimizer(B,S0,r,x0,OptAlg,lb,ub);
kappa=optimvars(1);
nubar=optimvars(2);
nu0=optimvars(3);
rho=optimvars(4);
chi=optimvars(5);
%}

%%%%%%%%%%%%%%%    PLOT AND PRINT RESULTS    %%%%%%%%%%%%%%%%%
%Plotter(kappa,nubar,nu0,rho,chi,S0,r,B,M,matur,SimPoints)
tic
Plotter_Sim(kappa,nubar,nu0,rho,chi,S0,r,M,T,repetitions,barr)

%Plotter3D(kappa,nubar,nu0,rho,chi,S0,r,B,b);
toc
tab=Printer(kappa,nubar,nu0,rho,chi,B,S0,r);
%openvar('tab')

beep


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%       FUNCTIONS       %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   PARAMETERS   %%%%
%kappa, nubar, nu0,rho,chi = Heston parameters
%S0=starting stock price
%r=risk-free rate
%T=maturity (in years)
%K=strike price
%L=number of simulation time steps
%M=number of simulated paths
%B=input data file
%x0=optimization starting parameters
%lb=parameters lower bounds
%ub=parameters upper bounds
%OptAlg=Optimization algorithm selected ("CMA" or "MultiStart")

%%%%%%%%%%%%%      DEFINE MINIMIZATION PROCEDURE      %%%%%%%%%%%%%%%
function optimvars=Optimizer(B,S0,r,x0,OptAlg,lb,ub)
fun=@(var)Hestoncal(var(1),var(2),var(3),var(4),var(5),S0,B,r,lb,ub); %function to be optimized.
%Hestoncal(kappa,nubar,nu0,rho,chi,S0,B,r,lb,ub)
tic
%MultiStart Algorithm
if OptAlg=="MultiStart"
    StartPoints=2; %Total number of starting points
    opts = optimoptions(@fmincon,'Display','off','Algorithm','active-set');
    problem = createOptimProblem('fmincon','objective',fun,'x0',x0,'lb',lb,'ub',ub,'options',opts);
    ms = MultiStart('UseParallel',true,'StartPointsToRun','bounds-ineqs','Display','off'); %Run minimizers in parallel
    [optimvars,f] = run(ms,problem,StartPoints);
    
    %CMA-ES Algorithm
elseif OptAlg=="CMA"
    optimvars=purecmaes(fun,x0);
    f=Hestoncal(optimvars(1),optimvars(2),optimvars(3),optimvars(4),optimvars(5),S0,B,r,lb,ub);  %optimal error value
end

%print optimization output
fprintf('kappa=%f, nubar=%f, nu0=%f\nrho=%f,   chi=%f\n%s, error=%f, t=%.0f\n\n',[optimvars(1),optimvars(2),optimvars(3),optimvars(4),optimvars(5),OptAlg,f,floor(toc/60)])
warning('on','MATLAB:integral:NonFiniteValue')
end



%%%  CALCULATE LEAST SQUARES ERROR BETWEEN MODEL AND DATA  IMPLIED VOL  %%%
function Total_Error=Hestoncal(kappa,nubar,nu0,rho,chi,S0,B,r,lb,ub)
%Check if bounds are being respected. If not, immediately output a big error so that test point is ignored
if kappa<=lb(1)||kappa>ub(1)||nubar<=lb(2)||nubar>ub(2)||nu0<=lb(3)||nu0>ub(3)||rho<=lb(4)||rho>ub(4)||chi<=lb(5)||chi>ub(5)
    Total_Error=1000;
    
    %If all boundaries are respected, calculate the error:
else
    LS=0;
    for i=1:size(B,1)
        %LS is the least squares error: (\sigma_market - \sigma_model)^2
        %Weight function: (1-|K/S_0-1|)^2
        LS=LS+(B(i,3)-HestonPrice(B(i,2),B(i,1),S0,r,kappa,nubar,nu0,rho,chi,"vol"))^2*(1-abs(B(i,2)-1))^2;
    end
    Total_Error=LS;
end
end



%OUTPUT THE PRICES/VOLATILITIES OF OPTIONS UNDER THE HESTON MODEL
%For more information check pages 6-9 of Cui et al.:  https://arxiv.org/pdf/1511.08718.pdf
function result=HestonPrice(K,T,S0,r,kappa,nubar,nu0,rho,chi,PriceVol)
warning('off','MATLAB:integral:NonFiniteValue')
%Calculate the price of options under Heston using the closed form solution
%We need the characteristic function, depicted in function "CharFuncHeston"
fun1=@(u)real(exp(-1i.*u.*log(K))./(1i.*u.*S0.*exp(r.*T)).*CharFuncHeston(u-1i,T,S0,r,kappa,nubar,nu0,rho,chi));
fun2=@(u)real(exp(-1i.*u.*log(K))./(1i.*u).*CharFuncHeston(u,T,S0,r,kappa,nubar,nu0,rho,chi));
P1=1/2+1/pi.*integral(fun1,0,400);
P2=1/2+1/pi.*integral(fun2,0,400);
if isnan(P1) || isinf(P1) || isnan(P2) || isinf(P2)
    result=10000;
else
    
    call=S0.*P1-exp(-r.*T).*K.*P2; %call option price
    
    
    
    if PriceVol=="price"   %if desired output is a price, return this result
        result=call;
    else                   %if desired output is a volatility, calculate implied volatility of price
        volatility=@(sigma)european_bs(S0,K,r,sigma,T,"call")-call;
        result=fzero(volatility,1);
    end
end
end



%OUTPUT THE CHARACTERISTIC FUNCTION FOR THE HESTON MODEL
%For more information check pages 6-9 of Cui et al.:  https://arxiv.org/pdf/1511.08718.pdf
function w=CharFuncHeston(u,T,S0,r,kappa,nubar,nu0,rho,chi)
xi=kappa-chi.*rho.*1i.*u;
d=sqrt(xi.^2+chi.^2.*(u.^2+1i.*u));
A1=(u.^2+1i.*u).*sinh(d.*T/2);
A2=d.*cosh(d.*T./2)+xi.*sinh(d.*T./2);
A=A1./A2;
D=log(d)+(kappa-d).*T./2-log((d+xi)./2+(d-xi)./2.*exp(-d.*T));
w=exp(1i.*u.*(log(S0)+r.*T)-T.*kappa.*nubar.*rho.*1i.*u./chi-nu0.*A+2*kappa.*nubar./chi.^2.*D);
end



function Plotter_Sim(kappa,nubar,nu0,rho,chi,S0,r,M,T,repetitions,barr)

K=0.4:0.01:1.6;

SimVol1=@(K,PriceVol)Pricer(kappa,nubar,nu0,rho,chi,K',S0,r,T,T*252*2,M,PriceVol,barr(1));
SimVol2=@(K,PriceVol)Pricer(kappa,nubar,nu0,rho,chi,K',S0,r,T,T*252*2,M,PriceVol,barr(2));
SimVol3=@(K,PriceVol)Pricer(kappa,nubar,nu0,rho,chi,K',S0,r,T,T*252*2,M,PriceVol,barr(3));
SimVolEuro=@(K,PriceVol)PricerEuro(kappa,nubar,nu0,rho,chi,K',S0,r,T,T*252*2,M,PriceVol);


    parfor j=1:repetitions
        DV_tmp1(j,:)=SimVol1(K,"vol");
        DV_tmp2(j,:)=SimVol2(K,"vol");
        DV_tmp3(j,:)=SimVol3(K,"vol");
        DVP_tmp1(j,:)=SimVol1(K,"price");
        DVP_tmp2(j,:)=SimVol2(K,"price");
        DVP_tmp3(j,:)=SimVol3(K,"price");
        DV_tmpEuro(j,:)=SimVolEuro(K,"vol");
        DVP_tmpEuro(j,:)=SimVolEuro(K,"price");
    end
    DV1=mean(DV_tmp1,1);
    DV2=mean(DV_tmp2,1);
    DV3=mean(DV_tmp3,1);
    DVP1=mean(DVP_tmp1,1);
    DVP2=mean(DVP_tmp2,1);
    DVP3=mean(DVP_tmp3,1);
    DVEuro=mean(DV_tmpEuro,1);
    DVPEuro=mean(DVP_tmpEuro,1);

    %Plot implied volatility function under the Heston model
    figure
    plot(K,DVEuro,'-.','LineWidth',2,'Color',[0.0010    0.60    0.8330]);
    hold on;
    plot(K,DV1,'LineWidth',2);
    plot(K,DV2,'--','LineWidth',2);
    plot(K,DV3,':','LineWidth',2);
    
    %Plot options
    xlim([0.4,1.6])
    ylim([0,1])
    box on;
    grid on;
    set(gca,'fontsize',12)
    xlabel('K/S_0');
    ylabel('\sigma_{imp} (yr^{-1/2})')
    pbaspect([1.5 1 1])


    lg={'European',['B=',num2str(barr(1))],['B=',num2str(barr(2))],['B=',num2str(barr(3))]};
    legend(lg,'Location','northeast','FontSize',11);

   
    figure
    plot(K,DVPEuro,'-.','LineWidth',2,'Color',[0.0010    0.60    0.8330]);
    hold on;
    plot(K,DVP1,'LineWidth',2);
    plot(K,DVP2,'--','LineWidth',2);
    plot(K,DVP3,':','LineWidth',2);
    
    %Plot options
    xlim([0.4,1.6])
    ylim([0,0.6])
    box on;
    grid on;
    set(gca,'fontsize',12)
    xlabel('K/S_0');
    ylabel('Option Price(€)')
    pbaspect([1.5 1 1])

    lg={'European Call',['B=',num2str(barr(1))],['B=',num2str(barr(2))],['B=',num2str(barr(3))]};
    legend(lg,'Location','northeast','FontSize',11);
end




%%%%%%%%%%%%%%    PLOT OPTIMIZATION RESULTS IN A SURFACE   %%%%%%%%%%%%%%%
function Plotter3D(kappa,nubar,nu0,rho,chi,S0,r,B)
figure

%Plot original data points
scatter3(B(:,2),B(:,1),B(:,3),30,'LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.3010    0.7450    0.9330]);
hold on;

%Plot implied volatility function under the Heston model
HestonVol=@(K,T)HestonPrice(K,T,S0,r,kappa,nubar,nu0,rho,chi,"vol");
[K,T] = meshgrid(0.4:0.01:1.6,0.5/12:0.1/12:0.5+0.5/12);
for i=1:size(K,1)
    for j=1:size(K,2)
        HV(i,j)=HestonVol(K(i,j),T(i,j));
    end
end
s=surf(K,T,HV);
%s.EdgeColor = 'interp';
s.EdgeAlpha=0.6;
s.FaceAlpha=0.85;
shading interp
hold on;

times=unique(B(:,1));   %array with all maturity dates
K2=0.4:0.01:1.6;
for i=1:size(times,1)
    for j=1:size(K2,2)
        HV2(j)=HestonVol(K2(j),times(i));
    end
    plot3(K2,ones(1,size(K2,2))*times(i),HV2,'LineWidth',2,'Color',[0.9500    0.200    0.1])
end

%Plot options
xlim([0.4,1.6])
ylim([0.5/12,0.5+0.5/12])
zlim([0 1])
box on;
grid on;
xlabel('K/S_0');
ylabel('T (days)');
zlabel('\sigma_{imp} (yr^{-1/2})')
yticks([1/12,2/12,3/12,4/12,5/12,6/12])
yticklabels({'21','42','63','84','105','126'})
%ylabel('\sigma_{imp} (yr^{-1/2})')
pbaspect([1 1.5 1])
%lgd=legend({'Market Data','Fitted Surface','Fitted Functions'},'FontSize',11);
%title(lgd,strcat(strcat("T=",num2str(T*252))," days"))
set(gca,'fontsize',11)
view(40,35)
M = view(gca);
R = M(1:3,1:3);
x = R*[1;0;0];
y = R*[0;1;0];
z = R*[0;0;1];
set(get(gca,'XLabel'),'rotation',360/(2*pi)*atan(x(2)/x(1)))
set(get(gca,'YLabel'),'rotation',360/(2*pi)*atan(y(2)/y(1)))



%Contour Plot
figure
contourf(K,T,HV,25)
pbaspect([1.5 1 1])
xlim([0.4,1.6])
ylim([0.5/12,0.5+0.5/12])
xlabel('K/S_0');
ylabel('T (days)');
yticks([1/12,2/12,3/12,4/12,5/12,6/12])
yticklabels({'21','42','63','84','105','126'})
box on;
colorbar;
set(gca,'fontsize',12)
end


%%% CALCULATE MONTE CARLO PRICE/IMPLIED VOLATILITY OF EUROPEAN OPTION UNDER HESTON %%%
%Options are assumed Call
%If output should be a price, PriceVol="price"
%If output should be an implied volatility, PriceVol="vol"
function Result=Pricer(kappa,nubar,nu0,rho,chi,C,S0,r,T,L,M,PriceVol,b)
dt = T/L;      %time steps
N=size(C,1);   %size of vector of volatilities to be output

S=S0*ones(M,1);      %define initial vector of stock prices
nu=nu0*ones(M,1);    %define the initial vector of volatilities
barrier = zeros(M,1);

for k = 1:L
    Z1=randn(M,1);                                 %vector of random variables
    Z2=rho*Z1+sqrt(1-rho^2)*randn(M,1);            %vector of random variables with correlation "rho" with vector Z1
    
    %Milstein discretization
    S(:)=S(:).*(1+r.*dt+sqrt(dt.*max(nu(:),0)).*Z1+0.5.*dt.*max(nu(:),0).*(Z1.^2-1));
    nu(:)=nu(:)+kappa.*(nubar-max(nu(:),0)).*dt+chi.*sqrt(dt.*max(nu(:),0)).*Z2+chi.^2./4.*dt.*(Z2.^2-1);  %Ensure variance never becomes negative
    barrier=max(barrier,S>b);
end
S=S.*barrier;


Y=zeros(M,N);    %matrix with paths' payoff  (for all inserted strikes)
for j=1:N
    for i=1:M
        Y(i,j) = max(S(i)-C(j,1),0);      %Calculate the payoff of all paths (assuming calls)
    end
end

if PriceVol=="price"
    Result=exp(-r*T)*mean(Y); %Output the discounted expected payoff
else
    Result=zeros(1,N);
    for j=1:N
        volatility=@(sigma)barrier_bs(S0,C(j,1),r,sigma,T,"call",b)-exp(-r*T)*mean(Y(:,j));
        res=fzero(volatility,0.01);   %Calculate the expected implied volatility
        
        if ~isnan(res)     %if no implied volatility is comaptible with the price, fzero outputs NaN
            Result(j)=res;
        else               %if no price is found, output zero implied vol
            Result(j)=0;
        end
    end
end

end



function Result=PricerEuro(kappa,nubar,nu0,rho,chi,C,S0,r,T,L,M,PriceVol)
dt = T/L;      %time steps
N=size(C,1);   %size of vector of volatilities to be output

S=S0*ones(M,1);      %define initial vector of stock prices
nu=nu0*ones(M,1);    %define the initial vector of volatilities

for k = 1:L
    Z1=randn(M,1);                                 %vector of random variables
    Z2=rho*Z1+sqrt(1-rho^2)*randn(M,1);            %vector of random variables with correlation "rho" with vector Z1
    
    %Milstein discretization
    S(:)=S(:).*(1+r.*dt+sqrt(dt.*max(nu(:),0)).*Z1+0.5.*dt.*max(nu(:),0).*(Z1.^2-1));
    nu(:)=nu(:)+kappa.*(nubar-max(nu(:),0)).*dt+chi.*sqrt(dt.*max(nu(:),0)).*Z2+chi.^2./4.*dt.*(Z2.^2-1);  %Ensure variance never becomes negative
end


Y=zeros(M,N);    %matrix with paths' payoff  (for all inserted strikes)
for j=1:N
    for i=1:M
        Y(i,j) = max(S(i)-C(j,1),0);      %Calculate the payoff of all paths (assuming calls)
    end
end

if PriceVol=="price"
    Result=exp(-r*T)*mean(Y); %Output the discounted expected payoff
else
    Result=zeros(1,N);
    for j=1:N
        volatility=@(sigma)european_bs(S0,C(j,1),r,sigma,T,"call")-exp(-r*T)*mean(Y(:,j));
        res=fzero(volatility,0.01);   %Calculate the expected implied volatility
        
        if ~isnan(res)     %if no implied volatility is comaptible with the price, fzero outputs NaN
            Result(j)=res;
        else               %if no price is found, output zero implied vol
            Result(j)=0;
        end
    end
end

end


%%% PRINT A TABLE WITH MODEL/MARKET IMPLIED VOL/PRICES AND REL. ERRORS %%%
function tab=Printer(kappa,nubar,nu0,rho,chi,B,S0,r)
format longG      %Change format for maximum precision
MKTVols=B(:,3);   %Market implied volatilities
MKTPrices=european_bs(S0,B(:,2),r,B(:,3),B(:,1),"call");  %Market (converted) prices

Vols=[]; Prices=[];
for i=1:size(B,1)
    vol=HestonPrice(B(i,2),B(i,1),S0,r,kappa,nubar,nu0,rho,chi,"vol");
    Vols=[Vols;vol];    %Model Implied volatilities
    Prices=[Prices;european_bs(S0,B(i,2),r,vol,B(i,1),"call")];  %Model Prices
end

%Output table
tab=[B(:,1)*252,B(:,2),MKTVols,Vols,abs(MKTVols-Vols)./MKTVols*100,MKTPrices,Prices,abs(MKTPrices-Prices)./MKTPrices*100];
format short
end



%%%%%%%%%%%%%%  CALCULATE BLACK-SCHOLES PRICE  %%%%%%%%%%%%%%%%%%%%
%If option is a call: putcall="call"
%If option is a put: putcall="put"
function price=european_bs(S0,K,r,sigma,T,putcall)
d1 = (log(S0./K) + (r + 0.5.*sigma.^2).*T)./(sigma.*sqrt(T));
d2 = d1 - sigma.*sqrt(T);
N1 = normcdf(d1);
N2 = normcdf(d2);
if putcall=="call"
    price = S0.*N1 - K.*exp(-r.*T).*N2;
elseif putcall=="put"
    price = S0.*N1 - K.*exp(-r*T).*N2 + K.*exp(-r.*T) - S0;
end
end

function price=barrier_bs(S0,K,r,sigma,T,putcall,B)
if B<K
    price=european_bs(S0,K,r,sigma,T,putcall);
else
x1=(B/S0)^(-1+2*r/sigma^2);
x2=(B/S0)^(1+2*r/sigma^2);
d3= (log(S0./B) + (r + 0.5.*sigma.^2).*T)./(sigma.*sqrt(T));
d4= (log(S0./B) - (r + 0.5.*sigma.^2).*T)./(sigma.*sqrt(T));
d5= (log(S0.*K./(B.^2)) - (r + 0.5.*sigma.^2).*T)./(sigma.*sqrt(T));
d6= (log(S0./B) + (r - 0.5.*sigma.^2).*T)./(sigma.*sqrt(T));
d7= (log(S0./B) - (r - 0.5.*sigma.^2).*T)./(sigma.*sqrt(T));
d8= (log(S0.*K./(B.^2)) - (r - 0.5.*sigma.^2).*T)./(sigma.*sqrt(T));
N3 = normcdf(d3);
N4 = normcdf(d4);
N5 = normcdf(d5);
N6 = normcdf(d6);
N7 = normcdf(d7);
N8 = normcdf(d8);
if putcall=="call"
    price = S0.*(N3+x2.*(N4-N5)) - K.*exp(-r.*T).*(N6+x1.*(N7-N8));
end
end
end


%%%%%%%%%%%%%%%      CMA-ES OPTIMIZATION ALGORITHM    %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%     FUNCTION CREATED BY HANSEN      %%%%%%%%%%%%%%%%%%%%
%%%%%%%%% ORIGINAL CODE: https://www.lri.fr/~hansen/purecmaes.m   %%%%%%%%%
function xmin=purecmaes(fun,x0)   % (mu/mu_w, lambda)-CMA-ES
% --------------------  Initialization --------------------------------
% User defined input parameters (need to be edited)
%  fun=@(var)SABRcal(var(1),var(2),var(3),beta,S0,Data,r);
%strfitnessfct = 'frosenbrock';  % name of objective/fitness function
N = size(x0,2);               % number of objective variables/problem dimension
%xmean = rand(N,1);    % objective variables initial point
xmean=x0';
sigma = 0.3;          % coordinate wise standard deviation (step size)
stopfitness = 1e-10;  % stop if fitness < stopfitness (minimization)
stopeval = 1e3*N^2;   % stop after stopeval number of function evaluations

% Strategy parameter setting: Selection
lambda = 4+floor(3*log(N));  % population size, offspring number
mu = lambda/2;               % number of parents/points for recombination
weights = log(mu+1/2)-log(1:mu)'; % muXone array for weighted recombination
mu = floor(mu);
weights = weights/sum(weights);     % normalize recombination weights array
mueff=sum(weights)^2/sum(weights.^2); % variance-effectiveness of sum w_i x_i

% Strategy parameter setting: Adaptation
cc = (4+mueff/N) / (N+4 + 2*mueff/N);  % time constant for cumulation for C
cs = (mueff+2) / (N+mueff+5);  % t-const for cumulation for sigma control
c1 = 2 / ((N+1.3)^2+mueff);    % learning rate for rank-one update of C
cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((N+2)^2+mueff));  % and for rank-mu update
damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; % damping for sigma
% usually close to 1
% Initialize dynamic (internal) strategy parameters and constants
pc = zeros(N,1); ps = zeros(N,1);   % evolution paths for C and sigma
B = eye(N,N);                       % B defines the coordinate system
D = ones(N,1);                      % diagonal D defines the scaling
C = B * diag(D.^2) * B';            % covariance matrix C
invsqrtC = B * diag(D.^-1) * B';    % C^-1/2
eigeneval = 0;                      % track update of B and D
chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));  % expectation of
%   ||N(0,I)|| == norm(randn(N,1))
% -------------------- Generation Loop --------------------------------
counteval = 0;  % the next 40 lines contain the 20 lines of interesting code
while counteval < stopeval
    
    % Generate and evaluate lambda offspring
    for k=1:lambda
        arx(:,k) = xmean + sigma * B * (D .* randn(N,1)); % m + sig * Normal(0,C)
        %arfitness(k) = feval(strfitnessfct, arx(:,k)); % objective function call
        arfitness(k) = fun(arx(:,k)); % objective function call
        counteval = counteval+1;
    end
    
    % Sort by fitness and compute weighted mean into xmean
    [arfitness, arindex] = sort(arfitness); % minimization
    xold = xmean;
    xmean = arx(:,arindex(1:mu))*weights;   % recombination, new mean value
    
    % Cumulation: Update evolution paths
    ps = (1-cs)*ps ...
        + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma;
    hsig = norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN < 1.4 + 2/(N+1);
    pc = (1-cc)*pc ...
        + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma;
    
    % Adapt covariance matrix C
    artmp = (1/sigma) * (arx(:,arindex(1:mu))-repmat(xold,1,mu));
    C = (1-c1-cmu) * C ...                  % regard old matrix
        + c1 * (pc*pc' ...                 % plus rank one update
        + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
        + cmu * artmp * diag(weights) * artmp'; % plus rank mu update
    
    % Adapt step size sigma
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
    
    % Decomposition of C into B*diag(D.^2)*B' (diagonalization)
    if counteval - eigeneval > lambda/(c1+cmu)/N/10  % to achieve O(N^2)
        eigeneval = counteval;
        C = triu(C) + triu(C,1)'; % enforce symmetry
        [B,D] = eig(C);           % eigen decomposition, B==normalized eigenvectors
        D = sqrt(diag(D));        % D is a vector of standard deviations now
        invsqrtC = B * diag(D.^-1) * B';
    end
    
    % Break, if fitness is good enough or condition exceeds 1e14, better termination methods are advisable
    if arfitness(1) <= stopfitness || max(D) > 1e7 * min(D)
        break;
    end
    
end % while, end generation loop

xmin = arx(:, arindex(1))'; % Return best point of last iteration.
% Notice that xmean is expected to be even
% better.
end