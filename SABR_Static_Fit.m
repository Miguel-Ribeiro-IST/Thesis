clear; %clear all variables. This prevents multiple simulations from interfering with one another
%Import data from file "Data_BNPP.txt" and assign it to matrix B. The file must be in the same folder
%1st column - maturities, 2nd column - strikes, 3rd column - implied volatilities
A = importdata('Data_BNPP.txt','\t',1);
B=A.data(:,:);

%%%%%%%%%%%%%%%%%%%%  INPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%
S0=17099.4;        %initial stock price
r = 0.06;          %risk-free rate. Forward prices in data file assumed r=0.06
matur=4;           %maturity at which we want to fit the data. If matur=5, the fifth maturity in the file is chosen.
OptMethod="CMA";
SimPoints=false;
ShowPrices=false;
x0=[0.5,-0.5,0.5,0.75];   %parameter starting values [alpha, rho, nu]
lb = [0,-1,0,0];       %parameter lower bounds
ub = [2,1,5,1];    %parameter upper bounds
M=50000;           %number of paths to be simulated
iterations=5;     %number of repetitions to be simulated (and then averaged)


%%%%%%%%%%%%%      ORIGINAL DATA MODIFICATIONS     %%%%%%%%%%%%
B(:,2)=B(:,2)/S0;     %normalize strike prices
S0=1;
B(:,1)=B(:,1)/252;    %convert maturities from days to years
times=unique(B(:,1));
T=times(matur);       %selected maturity
B=B(B(:,1)==T,:);     %only keep values of the maturity
L = T*252*2;          %number of steps in simulations



%%%%%%%%%%%%%%%%%%%%%    CALIBRATION      %%%%%%%%%%%%%%%%%%%%%%
%%{

optimvars=Optimizer(S0,B,r,x0,OptMethod,lb,ub);  %Obtain optimization variables
alpha=optimvars(1);
rho=optimvars(2);
nu=optimvars(3);
beta=optimvars(4);

%Plot optimization results
Plotter(alpha,rho,nu,beta,S0,r,T,L,M,iterations,B,OptMethod,SimPoints,ShowPrices)
beep
%}

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
function optimvars=Optimizer(S0,B,r,x0,MultiStoch,lb,ub)
fun=@(var)SABRcal(var(1),var(2),var(3),var(4),S0,B,r); %function to be optimized.
%alpha=var(1), rho=var(2), nu=var(3);
%SABRcal(alpha,rho,nu,beta,S0,B,r)

if MultiStoch=="MultiStart"
    rng default % For reproducibility
    opts = optimoptions(@fmincon,'Display','off','Algorithm','sqp');
    %problem = createOptimProblem('fmincon','objective',fun,'x0',x0,'lb',lb,'ub',ub,'options',opts,'nonlcon',@Feller_Condition);
    
    problem = createOptimProblem('fmincon','objective',fun,'x0',x0,'lb',lb,'ub',ub,'options',opts);
    
    ms = MultiStart('UseParallel',true,'StartPointsToRun','bounds-ineqs','Display','off');
    optimvars = run(ms,problem,10);
    %ms = GlobalSearch('StartPointsToRun','bounds-ineqs','Display','off');
    %[optimvars,f] = run(ms,problem);
    f=SABRcal(optimvars(1),optimvars(2),optimvars(3),optimvars(4),S0,B,r);
  
    
elseif MultiStoch=="CMA"
    optimvars=purecmaes(fun,x0);
    f=SABRcal(optimvars(1),optimvars(2),optimvars(3),optimvars(4),S0,B,r);
end

%print optimization output
fprintf(strcat("alpha=",strcat(num2str(optimvars(1)),strcat(",    rho=",strcat(num2str(optimvars(2)),strcat(",    nu=",strcat(num2str(optimvars(3)),strcat(",    beta=",strcat(num2str(optimvars(4)),strcat("\nerror=",strcat(num2str(f),"\n")))))))))));
end


%%%%%%%%%%%%%%%%%%%%%    PLOT OPTIMIZATION RESULTS    %%%%%%%%%%%%%%%%%%%%
function Plotter(alpha,rho,nu,beta,S0,r,T,L,M,iterations,B,OptMethod,SimPoints,ShowPrices)
figure
SABRVol=@(K)sigmaSABR(alpha,rho,nu,beta,K,S0*exp(r*T),T);   %SABR implied volatility as a function of K
if SimPoints && ShowPrices
    ax1 = subplot(1,2,1);
else ax1 = subplot(1,1,1);
end

scatter(ax1,B(:,2),B(:,3),'.');
hold on;
fplot(ax1,SABRVol,[min(B(:,2)) max(B(:,2))],'k')


if SimPoints
    for i=1:size(B,1)
        if ShowPrices
            P(i)=european_bs(S0,B(i,2),r,B(i,3),B(i,1),"call");  %obtain the BS price from each of the data's implied volatilities
            %european_bs(S0,K,r,sigma0,T,putcall)
            %Obtain Monte Carlo prices and implied volatilities
            Price(i)=Pricer(alpha,rho,nu,beta,B(i,2),S0,r,T,L,M,iterations,"price");
        end
        
        Volatility(i)=Pricer(alpha,rho,nu,beta,B(i,2),S0,r,T,L,M,iterations,"vol");
        %Pricer(alpha,rho,nu,beta,K,f,r,T,L,M,iterations,PriceVol)
        %Plot prices from data and from Monte-Carlo
    end
    hold on;
    scatter(ax1,B(:,2),Volatility(:),'x');
    
    if ShowPrices
        ax2 = subplot(1,2,2);
        scatter(ax2,B(:,2),P(:),'.');
        hold on;
        scatter(ax2,B(:,2),Price(:),'x');
        
        %Show fit results in plot titles
        text=strcat(strcat(strcat(strcat(strcat("\beta=",num2str(beta)),strcat(", paths=",num2str(M))),strcat(", iterations=",num2str(iterations))),strcat(", maturity=",num2str(T*252))),strcat(", method=",OptMethod));
        text2=strcat(strcat(strcat("\alpha=",num2str(alpha)),strcat(", \rho=",num2str(rho))),strcat(", \nu=",num2str(nu)));
        title(ax1,{text,'Volatilities'})
        title(ax2,{text2,'Prices'})
    end
else
    text=strcat(strcat(strcat(strcat(strcat(strcat("\beta=",num2str(beta)),strcat(", paths=",num2str(M))),strcat(", iterations=",num2str(iterations))),strcat(", maturity=",num2str(T*252))),strcat(", method=",OptMethod)),strcat(strcat(strcat(", \alpha=",num2str(alpha)),strcat(", \rho=",num2str(rho))),strcat(", \nu=",num2str(nu))));
    title(ax1,{text,'Volatilities'})
end

%Plot implied volatilities from data, from Monte-Carlo and the SABR implied volatility function




end


%%%  CALCULATE LEAST SQUARES ERROR BETWEEN MODEL AND DATA  IMPLIED VOL  %%%
function Total_Error=SABRcal(alpha,rho,nu,beta,S0,B,r)
if alpha<0 || alpha>5 || rho<-1 || rho>1 || nu<0 || nu>5 || beta<0 || beta>1
    Total_Error=10000;
else
LS=0;
for i=1:size(B,1)
    %LS is the least squares error
    %Function sigmaSABR outputs the error between model and data implied
    %volatilities using the original Hagan formula
    LS=LS+(((B(i,3)-sigmaSABR(alpha,rho,nu,beta,B(i,2),S0*exp(r*B(i,1)),B(i,1)))./B(i,3))^2);
    %sigmaSABR(alpha,rho,nu,beta,K,f,T)
end
Total_Error=LS;
end
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
function Result_Avg=Pricer(alpha,rho,nu,beta,K,S0,r,T,L,M,iterations,PriceVol)
dt = T/L;      %time steps

%%NOTE:
%%%%We don't need to simulate an entire matrix of forwards for all time steps and for all paths.
%%%%We just need to simulate one vector of forwards for all paths and update it at each time step
parfor iter=1:iterations    %Perform the "for" cycle in parallel
    S = S0*ones(M,1);        %define initial vector of forwards
    sigma=alpha*ones(M,1);    %define the initial vector of volatilities
    
    for k = 1:L
        Z1=randn(M,1);                                 %vector of random variables
        Z2=rho*Z1+sqrt(1-rho^2)*randn(M,1);            %vector of random variables with correlation "rho" with vector Z1
        %v=alp.*(F.^(beta-1));                          %intermediate variable
        %F(:)=F(:).*exp(v.*Z2*sqrt(dt)-v.^2*dt/2);      %update forwards vector
        %alp(:)=alp(:).*exp(nu*sqrt(dt)*Z1-nu^2*dt/2);  %update volatilities vector
        S(:)=max(S(:).*(1+r*dt)+exp(-r*(T-dt*k)*(1-beta)).*sigma(:).*S(:).^beta.*sqrt(dt).*Z1+beta/2*exp(-2*r*(T-dt*k)*(1-beta))*sigma(:).^2.*S(:).^(2*beta-1)*dt.*(Z1.^2-1),0);
        sigma(:)=max(sigma(:).*(1+nu*sqrt(dt).*Z2+nu^2/2*dt*(Z2.^2-1)),0);
    end
    
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