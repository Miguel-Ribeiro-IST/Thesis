clearvars -except kappa nubar nu0 rho chi; %clear all variables. This prevents multiple simulations from interfering with one another
%Import data from file "Data_BNPP_2.txt" and assign it to matrix B. The file must be in the same folder
%1st column - maturities, 2nd column - strikes, 3rd column - implied volatilities
A = importdata('Data_BNPP.txt','\t',1);
B=A.data(:,:);

%%%%%%%%%%%%%%%%%%%%  INPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%
S0=17099.4;        %initial stock price
r = 0;             %risk-free rate. Forward prices in data file assumed r=0.0
matur=4;           %maturity until which we want to fit the data. If matur=5, all maturities until the fifth are chosen.
OptAlg="MultiStart";      %"CMA" or "MultiStart" optimization algorithms

%%%%%%%%%%%%%%%%%%%   MONTE CARLO SIMULATION %%%%%%%%%%%%%%
%After calibrating all the model's parameters, we may simulate the implied volatilities using Monte Carlo
SimPoints=false;   %true or false - define if Monte Carlo simulation should be executed
M=50000;           %number of paths to be simulated
iterations=5;      %number of repetitions to be simulated (and then averaged) to reduce error


%%%%%%%%%%%%%      ORIGINAL DATA MODIFICATIONS     %%%%%%%%%%%%
B(:,2)=B(:,2)/S0;                %normalize strike prices
S0=1;
B(:,1)=B(:,1)/252;               %convert maturities from days to years
times=unique(B(:,1));
B=B(B(:,1)<=times(matur),:);     %only keep values of the maturity

P=B;  %Generate matrix with European call prices instead of implied vol.
for i=1:size(B,1)
    P(i,3)=european_bs(S0,B(i,2),r,B(i,3),B(i,1),"call");
end


%%%%%%%%%%%%%%%%%%%%%    CALIBRATION      %%%%%%%%%%%%%%%%%%%%%%
%    [kappa, nubar, nu0, rho, chi]
x0=  [0.01,  1, 1, 0.9, 0.001];    %parameter starting values
lb = [0,     0, 0, -1,  0];          %parameter lower bounds
ub = [50,    2, 2, 1,   5];          %parameter upper bounds
%%{
optimvars=Optimizer(B,S0,r,x0,OptAlg,lb,ub);
kappa=optimvars(1);
nubar=optimvars(2);
nu0=optimvars(3);
rho=optimvars(4);
chi=optimvars(5);
%}

%%%%%%%%%%%%%%%%%%%    PLOT RESULTS    %%%%%%%%%%%%%%%%%%%%%
Plotter(kappa,nubar,nu0,rho,chi,S0,r,B,M,iterations,matur,SimPoints)

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
%M=number of simulated paths
%iterations=number of repetitions of the simulations (to reduce error)
%B=input data file
%x0=optimization starting parameters
%lb=parameters lower bounds
%ub=parameters upper bounds
%OptAlg=Optimization algorithm selected ("CMA" or "MultiStart")


%OUTPUT THE PRICES/VOLATILITIES OF OPTIONS UNDER THE HESTON MODEL
%For more information check pages 6-9 of Cui et al.:  https://arxiv.org/pdf/1511.08718.pdf
function result=HestonPrice(K,T,S0,r,kappa,nubar,nu0,rho,chi,PriceVol)
%Calculate the price of options under Heston using the closed form solution
%We need the characteristic function, depicted in function "CharFuncHeston"
fun1=@(u)real(exp(-1i.*u.*log(K))./(1i.*u.*S0.*exp(r.*T)).*CharFuncHeston(u-1i,T,S0,r,kappa,nubar,nu0,rho,chi));
fun2=@(u)real(exp(-1i.*u.*log(K))./(1i.*u).*CharFuncHeston(u,T,S0,r,kappa,nubar,nu0,rho,chi));
P1=1/2+1/pi.*integral(fun1,0,300);
P2=1/2+1/pi.*integral(fun2,0,300);
call=S0.*P1-exp(-r.*T).*K.*P2; %call option price

if PriceVol=="price"   %if desired output is a price, return this result
    result=call;
else                   %if desired output is a volatility, calculate implied volatility of price
    volatility=@(sigma)european_bs(S0,K,r,sigma,T,"call")-call;
    result=fzero(volatility,0.25);
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


%%%%%%%%%%%%%      DEFINE MINIMIZATION PROCEDURE      %%%%%%%%%%%%%%%
function optimvars=Optimizer(B,S0,r,x0,OptAlg,lb,ub)
fun=@(var)Hestoncal(var(1),var(2),var(3),var(4),var(5),S0,B,r,lb,ub); %function to be optimized.
%Hestoncal(kappa,nubar,nu0,rho,chi,S0,P,r)
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
    f=Hestoncal(optimvars(1),optimvars(2),optimvars(3),optimvars(4),optimvars(5),S0,B,r,lb,ub);
end

%print optimization output
fprintf('kappa=%f, nubar=%f, nu0=%f\nrho=%f,   chi=%f\nOptAlg=%s, error=%f, t=%.0f\n\n',[optimvars(1),optimvars(2),optimvars(3),optimvars(4),optimvars(5),OptAlg,f,floor(toc/60)])
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




%%%%%%%%%%%%%%%%%%%%%    PLOT OPTIMIZATION RESULTS    %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%    PLOT MULTIPLE MATURITIES IN A SINGLE FIGURE  %%%%%%%%%%%%
function MultiPlotter(kappa,nubar,nu0,rho,chi,S0,r,B,M,iterations,matur,SimPoints,OptAlg)
figure
times=unique(B(:,1));   %array with all maturity dates

for iter=1:matur
    if matur>1
    ax(iter) = subplot(2,ceil(matur/2),iter);  %divide figure into two rows
    else
    ax(iter) = subplot(1,1,1);
    end
    T=times(iter);
    
    C=B(B(:,1)==T,2:3);  %Implied volatilities for the selected maturity
    
    
    
    %If the user chose to use Monte Carlo, after calibration, to check model validity, calculate implied volatilities under MC
    if SimPoints
        for i=1:size(C,1)
            Volatility(i)= Pricer(kappa,nubar,nu0,rho,chi,C(i,1),S0,r,T,T*252*2,M,iterations,"vol");
        end
        scatter(ax(iter),C(:,1),Volatility(:),'x');
        hold on;
    end
    
    %Plot original data points
    scatter(ax(iter),C(:,1),C(:,2),'.');
    hold on;
    
    %Plot implied volatility function under the Heston model
    HestonVol=@(K)HestonPrice(K,T,S0,r,kappa,nubar,nu0,rho,chi,"vol");
    fplot(ax(iter),HestonVol,[0.4,1.6])
    hold on;
    
    %Plot options
    xlim([0.4,1.6])
    ylim([0.2,1])
    box on;
    grid on;
    title(ax(iter),strcat(strcat(strcat(num2str(T*252)," days  ("),num2str(T*252/21))," months)"))
    
    clear Volatility
end

%Show parameter values in plots as titles
if matur>1
text1=strcat(strcat(strcat(num2str(times(1)*252)," days  ("),num2str(times(1)*252/21))," months)");
if SimPoints
    vars1=strcat(strcat(strcat(strcat("paths=",num2str(M)),strcat(",  iterations=",num2str(iterations))),",  method="),OptAlg);
else
    vars1=strcat("method=",OptAlg);
end
text2=strcat(strcat(strcat(num2str(times(2)*252)," days  ("),num2str(times(2)*252/21))," months)");
vars2=strcat("\kappa=",strcat(num2str(kappa),strcat(", \nu_{mean}=",strcat(num2str(nubar),strcat(", \nu_0=",strcat(num2str(nu0),strcat(strcat(strcat(", \rho=",num2str(rho)),(strcat(", \chi=",num2str(chi)))))))))));
title(ax(1),{vars1,text1})
title(ax(2),{vars2,text2})
end

end



%%%%%%%%%%%%%%%%%%%%%    PLOT OPTIMIZATION RESULTS    %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%    PLOT MATURITIES IN DIFFERENT FIGURES  %%%%%%%%%%%%
function Plotter(kappa,nubar,nu0,rho,chi,S0,r,B,M,iterations,matur,SimPoints)

times=unique(B(:,1));   %array with all maturity dates

for iter=1:matur
    figure
    T=times(iter);
    
    C=B(B(:,1)==T,2:3);  %Implied volatilities for the selected maturity
    
    
    
    %If the user chose to use Monte Carlo, after calibration, to check model validity, calculate implied volatilities under MC
    if SimPoints
        for i=1:size(C,1)
            Volatility(i)= Pricer(kappa,nubar,nu0,rho,chi,C(i,1),S0,r,T,T*252*2,M,iterations,"vol");
        end
        scatter(C(:,1),Volatility(:),100,'x','LineWidth',1.5);
        hold on;
    end
    
    %Plot original data points
    scatter(C(:,1),C(:,2),100,'x','LineWidth',1.5);
    hold on;
    
    %Plot implied volatility function under the Heston model
    HestonVol=@(K)HestonPrice(K,T,S0,r,kappa,nubar,nu0,rho,chi,"vol");
    fplot(HestonVol,[0.4,1.6],'LineWidth',1.5)
    hold on;
    
    %Plot options
    xlim([0.4,1.6])
    ylim([0.2,1])
    box on;
    grid on;
    set(gca,'fontsize',12)
    xlabel('K/S_0');
    ylabel('\sigma_{imp} (yr^{-1/2})')
    pbaspect([1.5 1 1])
    lgd=legend({'Market Data','Fitted Function'},'Location','northeast','FontSize',11);
    title(lgd,strcat(strcat("T=",num2str(T*252))," days"))
    
    clear Volatility
end

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
    %x = zeros(M,1);        %define initial vector of forwards
    S=S0*ones(M,1);
    nu=nu0*ones(M,1);    %define the initial vector of volatilities
    
    for k = 1:L
        Z1=randn(M,1);                                 %vector of random variables
        Z2=rho*Z1+sqrt(1-rho^2)*randn(M,1);            %vector of random variables with correlation "rho" with vector Z1
        
        %Milstein discretization
        S(:)=S(:).*(1+r.*dt+sqrt(dt.*nu(:)).*Z1+0.5.*dt.*nu(:).*(Z1.^2-1));
        nu(:)=max(nu(:)+kappa.*(nubar-nu(:)).*dt+chi.*sqrt(dt.*nu(:)).*Z2+chi.^2./4.*dt.*(Z2.^2-1),0);  %Ensure variance never becomes negative
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