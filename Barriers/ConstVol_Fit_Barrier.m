clear;
A = importdata('Data_BNPP.txt','\t',1);
B=A.data(:,:);




%%%%%%%%%%%%%%%% CHECK MATURITY DATE AND RESPECTIVE PARAMETERS %%%%%%%%%





%matur=1;sigma=0.3174;
matur=2;sigma=0.2918;
%matur=3;sigma=0.2742;
%matur=4;sigma=0.2518;


%%%%%%%%%%%%%%%%%%%%  INPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%
S0=17099.4;
r = 0;
OptAlg="CMA";

%%%%%%%%%%%%%%%%%%%   MONTE CARLO SIMULATION %%%%%%%%%%%%%%
SimPoints=false;
M=10000;
repetitions=10;
barr=[1.05,1.1,1.2];
%L=T*252*2


%%%%%%%%%%%%%      ORIGINAL DATA MODIFICATIONS     %%%%%%%%%%%%
B(:,2)=B(:,2)/S0;
S0=1;
B(:,1)=B(:,1)/252;
times=unique(B(:,1));
T=times(matur);
B=B(B(:,1)==T,:);


%%%%%%%%%%%%%%%%%%%%%    CALIBRATION      %%%%%%%%%%%%%%%%%%%%%%
Plotter_Sim(sigma,S0,r,M,T,repetitions,barr)


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
%iterations=number of repetitions of the simulations (to reduce error)
%B=input data file
%x0=optimization starting parameters
function Plotter_Sim(sigma,S0,r,M,T,repetitions,barr)

K=0.4:0.01:1.6;

SimVol1=@(K,PriceVol)Pricer(sigma,K',S0,r,T,T*252*2,M,PriceVol,barr(1));
SimVol2=@(K,PriceVol)Pricer(sigma,K',S0,r,T,T*252*2,M,PriceVol,barr(2));
SimVol3=@(K,PriceVol)Pricer(sigma,K',S0,r,T,T*252*2,M,PriceVol,barr(3));
SimVolEuro=@(K,PriceVol)PricerEuro(sigma,K',S0,r,T,T*252*2,M,PriceVol);



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
    lgd=legend(lg,'Location','northeast','FontSize',11);
    title(lgd,"Constant Vol")
   
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

    lg={'European',['B=',num2str(barr(1))],['B=',num2str(barr(2))],['B=',num2str(barr(3))]};
    lgd=legend(lg,'Location','northeast','FontSize',11);
    title(lgd,"Constant Vol")
end


function Result=Pricer(sigma,C,S0,r,T,L,M,PriceVol,b)
dt = T/L;
N=size(C,1);

S = S0*ones(M,1);
barrier = zeros(M,1);
for k = 1:L
    %Euler-Maruyama discretization
    S(:)=S(:).*(1+r*dt+sigma*sqrt(dt)*randn(M,1));
    barrier=max(barrier,S>b);
end
S=S.*barrier;

Y=zeros(M,N);

for j=1:N
    for i=1:M
        Y(i,j) = max(S(i)-C(j,1),0);
    end
end

if PriceVol=="price"
    Result=exp(-r*T)*mean(Y);
else
    Result=zeros(1,N);
    for j=1:N
        volatility=@(sigma)barrier_bs(S0,C(j,1),r,sigma,T,"call",b)-exp(-r*T)*mean(Y(:,j));
        res=fzero(volatility,0.0001);
        
        if ~isnan(res)
            Result(j)=res;
        else
            Result(j)=0;
        end
    end
end
end


function Result=PricerEuro(sigma,C,S0,r,T,L,M,PriceVol)
dt = T/L;
N=size(C,1);

S = S0*ones(M,1);
for k = 1:L
    %Euler-Maruyama discretization
    S(:)=S(:).*(1+r*dt+sigma*sqrt(dt)*randn(M,1));
end

Y=zeros(M,N);

for j=1:N
    for i=1:M
        Y(i,j) = max(S(i)-C(j,1),0);
    end
end

if PriceVol=="price"
    Result=exp(-r*T)*mean(Y);
else
    Result=zeros(1,N);
    for j=1:N
        volatility=@(sigma)european_bs(S0,C(j,1),r,sigma,T,"call")-exp(-r*T)*mean(Y(:,j));
        res=fzero(volatility,0.0001);
        
        if ~isnan(res)
            Result(j)=res;
        else
            Result(j)=0;
        end
    end
end
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