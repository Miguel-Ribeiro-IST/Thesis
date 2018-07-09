clearvars -except alpha rho nu beta; %clear all variables. This prevents multiple simulations from interfering with one another
%Import data from file "Data_BNPP.txt" and assign it to matrix B. The file must be in the same folder
%1st column - maturities, 2nd column - strikes, 3rd column - implied volatilities
A = importdata('Data_BNPP.txt','\t',1);
B=A.data(:,:);

%%%%%%%%%%%%%%%%%%%%  INPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%
S0=17099.4;               %initial stock price
r = 0;                    %risk-free rate. Forward prices in data file assumed r=0
matur=4;                  %maturity at which we want to fit the data. If matur=5, only the fifth maturity in the file is chosen.
OptAlg="CMA";      %"CMA" or "MultiStart" optimization algorithms


%%%%%%%%%%%%%%%%%%%   MONTE CARLO SIMULATION %%%%%%%%%%%%%%
%After calibrating all the model's parameters, we may want to simulate the implied volatilities using Monte Carlo
SimPoints=false;   %true or false - define if Monte Carlo simulation should be executed
M=100000;           %number of paths to be simulated
repetitions=100;
%L=T*252*2


%%%%%%%%%%%%%      ORIGINAL DATA MODIFICATIONS     %%%%%%%%%%%%
B(:,2)=B(:,2)/S0;     %normalize strike prices
S0=1;
B(:,1)=B(:,1)/252;    %convert maturities from days to years
times=unique(B(:,1));
T=times(matur);       %selected maturity
B=B(B(:,1)==T,:);     %only keep values of the maturity


%%%%%%%%%%%%%%%%%%%%%    CALIBRATION      %%%%%%%%%%%%%%%%%%%%%%
%    [alpha, rho, nu, beta]
x0 = [0.5, -0.5, 0.5, 0.75];   %parameter starting values
lb = [0,   -1,   0,   0];       %parameter lower bounds
ub = [2,    1,   5,   1];    %parameter upper bounds
%%{
optimvars=Optimizer(S0,B,r,x0,OptAlg,lb,ub);  %Obtain optimization variables
alpha=optimvars(1);
rho=optimvars(2);
nu=optimvars(3);
beta=optimvars(4);
%}
tic
%%%%%%%%%%%%%%%%%%%    PLOT RESULTS    %%%%%%%%%%%%%%%%%%%%%
%Plotter(alpha,rho,nu,beta,S0,r,T,M,B,SimPoints)
Plotter_Sim(alpha,rho,nu,beta,S0,r,T,M,B,repetitions)
toc
tab=Printer(alpha,rho,nu,beta,B,S0,r,T);
%openvar('tab')

beep

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
%B=input data file
%x0=optimization starting parameters
%lb=parameters lower bounds
%ub=parameters upper bounds
%OptAlg=Optimization algorithm selected ("CMA" or "MultiStart")

%%%%%%%%%%%%%      DEFINE MINIMIZATION PROCEDURE      %%%%%%%%%%%%%%%
function optimvars=Optimizer(S0,B,r,x0,OptAlg,lb,ub)
fun=@(var)SABRcal(var(1),var(2),var(3),var(4),S0,B,r,lb,ub); %function to be optimized.
%SABRcal(alpha,rho,nu,beta,S0,B,r,lb,ub)

%MultiStart Algorithm
if OptAlg=="MultiStart"
    StartPoints=2; %Total number of starting points
    opts = optimoptions(@fmincon,'Display','off','Algorithm','sqp');
    problem = createOptimProblem('fmincon','objective',fun,'x0',x0,'lb',lb,'ub',ub,'options',opts);
    ms = MultiStart('UseParallel',true,'StartPointsToRun','bounds-ineqs','Display','off'); %Run minimizers in parallel
    [optimvars,f] = run(ms,problem,StartPoints);
    
    
    %CMA-ES Algorithm
elseif OptAlg=="CMA"
    optimvars=purecmaes(fun,x0);
    f=SABRcal(optimvars(1),optimvars(2),optimvars(3),optimvars(4),S0,B,r,lb,ub);  %optimal error value
end


%print optimization output
fprintf('alpha=%f, rho=%f\nnu=%f,    beta=%f\n%s, error=%f\n\n',[optimvars(1),optimvars(2),optimvars(3),optimvars(4),OptAlg,f])
end


%%%  CALCULATE LEAST SQUARES ERROR BETWEEN MODEL AND DATA  IMPLIED VOL  %%%
function Total_Error=SABRcal(alpha,rho,nu,beta,S0,B,r,lb,ub)
%Check if bounds are being respected. If not, immediately output a big error so that test point is ignored
if alpha<lb(1) || alpha>ub(1) || rho<lb(2) || rho>ub(2) || nu<lb(3) || nu>ub(3) || beta<lb(4) || beta>ub(4)
    Total_Error=1000;
    
    %If all boundaries are respected, calculate the error:
else
    LS=0;
    for i=1:size(B,1)
        %LS is the least squares error: (\sigma_market - \sigma_model)^2
        %Weight function: (1-|K/S_0-1|)^2
        LS=LS+(B(i,3)-sigmaSABR(alpha,rho,nu,beta,B(i,2),S0*exp(r*B(i,1)),B(i,1)))^2*(1-abs(B(i,2)-1))^2;
    end
    Total_Error=LS;
end
end


%%%%%%%%%%% CALCULATE SABR MODEL IMPLIED VOLATILITY FROM HAGAN  %%%%%%%%%%
%%Find function in page 9 of file github.com/Miguel-Ribeiro-IST/Thesis/blob/master/References/Hagan_SABR.pdf
function sigma=sigmaSABR(alpha,rho,nu,beta,K,f,T)
z=nu./alpha.*(f.*K).^((1-beta)./2).*log(f./K);
x=log((sqrt(1-2.*rho.*z+z.^2)+z-rho)./(1-rho));
z(z==0,:)=1;
x(x==0,:)=1;

if beta==0
    sigma=alpha.*log(f./K)./(f-K).*(z./x).*(1+T.*(alpha.^2./(24.*f.*K)+(2-3*rho.^2)./24.*nu.^2));
elseif beta==1
    sigma=alpha.*(z./x).*(1+T.*(1/4*rho.*alpha.*nu+1/24*(2-3*rho.^2).*nu.^2));
else
    sigma=alpha./((f.*K).^((1-beta)./2).*(1+(1-beta).^2/24.*(log(f./K)).^2+(1-beta).^4./1920.*(log(f./K)).^4)).*(z./x).*(1+((1-beta).^2./24.*alpha.^2./(f.*K).^(1-beta)+1/4.*rho.*nu.*beta.*alpha./(f.*K).^((1-beta)./2)+(2-3.*rho.^2)./24.*nu.^2).*T);
end
end


%%%%%%%%%%%%%%%%%%%%%    PLOT OPTIMIZATION RESULTS    %%%%%%%%%%%%%%%%%%%%
function Plotter(alpha,rho,nu,beta,S0,r,T,M,B,SimPoints)
figure

%If the user chose to use Monte Carlo, after calibration, to check model validity, calculate implied volatilities under MC
if SimPoints
    Volatility=Pricer(alpha,rho,nu,beta,B(:,2:3),S0,r,T,T*252*2,M,"vol");
    scatter(B(:,2),Volatility(:),100,'+','LineWidth',1.5);
    hold on;
end

%Plot original data points
scatter(B(:,2),B(:,3),100,'x','LineWidth',1.5);
hold on;

%Plot implied volatility function under the SABR model
SABRVol=@(K)K.*0+sigmaSABR(alpha,rho,nu,beta,K',S0.*exp(r.*T),T)';
fplot(SABRVol,[0.4,1.6],'LineWidth',1.5)

%Plot options
xlim([0.4,1.6])
ylim([0.2,1])
box on;
grid on;
set(gca,'fontsize',12)
xlabel('K/S_0');
ylabel('\sigma_{imp} (yr^{-1/2})')
pbaspect([1.5 1 1])
if SimPoints
    lgd=legend({'Simulated Volatilities','Market Data','Fitted Function'},'Location','northeast','FontSize',11);
else
    lgd=legend({'Market Data','Fitted Function'},'Location','northeast','FontSize',11);
end
title(lgd,strcat(strcat("T=",num2str(T*252))," days"))
end


function Plotter_Sim(alpha,rho,nu,beta,S0,r,T,M,B,repetitions)
figure

%If the user chose to use Monte Carlo, after calibration, to check model validity, calculate implied volatilities under MC
K=(0.4:0.01:1.6);
SimVol=@(K)Pricer(alpha,rho,nu,beta,K,S0,r,T,T*252*2,M,"vol");
for j=1:repetitions
    Mdl_tmp(j,:)=SimVol(K');
end
Mdl=mean(Mdl_tmp);
Mdlmax90=quantile(Mdl_tmp,0.95,1);
Mdlmin10=quantile(Mdl_tmp,0.05,1);


%Plot original data points
scatter(B(:,2),B(:,3),100,[0    0.1470    0.6410],'x','LineWidth',1.5);
hold on;

%Plot implied volatility function under the SABR model
SABRVol=@(K)K.*0+sigmaSABR(alpha,rho,nu,beta,K',S0.*exp(r.*T),T)';
fplot(SABRVol,[0.4,1.6],'LineWidth',2,'Color',[1.00    0.3050    0.0580])
hold on;

plot(K,Mdl,'-.','LineWidth',2,'Color',[0.0010    0.60    0.8330]);


hold on;
K2 = [K, fliplr(K)]; % Use ; instead of ,
inBetween = [Mdlmax90, fliplr(Mdlmin10)]; % Use ; instead of ,
fill(K2, inBetween,[0    0.150    0.830],'FaceAlpha',0.2,'EdgeAlpha',0);

%Plot options
xlim([0.4,1.6])
ylim([0,1])
box on;
grid on;
set(gca,'fontsize',12)
xlabel('K/S_0');
ylabel('\sigma_{imp} (yr^{-1/2})')
pbaspect([1.5 1 1])


h = get(gca,'Children');
lgd=legend([h(4) h(3) h(2) h(1)],{'Market Data','Theoretical Function','Simulated Function (mean)','95% Confidence Interval'},'Location','northeast','FontSize',11);
title(lgd,strcat(strcat("T=",num2str(T*252))," days"))
set(gca,'Children',[h(4) h(2) h(3) h(1)])
end


%%% CALCULATE MONTE CARLO PRICE/IMPLIED VOLATILITY OF EUROPEAN OPTION UNDER SABR %%%
%Options are assumed Call
%If output should be a price, PriceVol="price"
%If output should be an implied volatility, PriceVol="vol"
function Result=Pricer(alpha,rho,nu,beta,C,S0,r,T,L,M,PriceVol)
dt = T/L;      %time steps
N=size(C,1);   %size of vector of volatilities to be output


S = S0*ones(M,1);        %define initial vector of forwards
sigma=alpha*ones(M,1);    %define the initial vector of volatilities

for k = 1:L
    Z1=randn(M,1);                                 %vector of random variables
    Z2=rho*Z1+sqrt(1-rho^2)*randn(M,1);            %vector of random variables with correlation "rho" with vector Z1
    
    %Milstein discretization
    S(:)=max(S(:),0).*(1+r*dt)+exp(-r*(T-dt*k)*(1-beta)).*sigma(:).*max(S(:),0).^beta.*sqrt(dt).*Z1+beta/2*exp(-2*r*(T-dt*k)*(1-beta))*sigma(:).^2.*max(S(:),0).^(2*beta-1)*dt.*(Z1.^2-1);
    sigma(:)=sigma(:).*(1+nu*sqrt(dt).*Z2+nu^2/2*dt*(Z2.^2-1));
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
        res=fzero(volatility,0.0001);   %Calculate the expected implied volatility
        
        if ~isnan(res)     %if no implied volatility is comaptible with the price, fzero outputs NaN
            Result(j)=res;
        else               %if no price is found, return 0
            Result(j)=0;
        end
    end
end

end



function tab=Printer(alpha,rho,nu,beta,B,S0,r,T)
format longG
MKTVols=B(:,3);
MKTPrices=european_bs(S0,B(:,2),r,B(:,3),B(:,1),"call");

SABRVol=@(K)sigmaSABR(alpha,rho,nu,beta,K,S0.*exp(r.*T),T);
Vols=SABRVol(B(:,2));
Prices=european_bs(S0,B(:,2),r,Vols,B(:,1),"call");


tab=[B(:,2),MKTVols,Vols,abs(MKTVols-Vols)./MKTVols*100,MKTPrices,Prices,abs(MKTPrices-Prices)./MKTPrices*100];
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