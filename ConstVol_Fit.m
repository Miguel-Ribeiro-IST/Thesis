clear; %clear all variables. This prevents multiple simulations from interfering with one another
%Import data from file "Data_BNPP.txt" and assign it to matrix B. The file must be in the same folder
%1st column - maturities, 2nd column - strikes, 3rd column - implied volatilities
A = importdata('Data_BNPP_2.txt','\t',1);
B=A.data(:,:);

%%%%%%%%%%%%%%%%%%%%  INPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%
S0=17099.4;        %initial stock price
r = 0;          %risk-free rate. Forward prices in data file assumed r=0.06
matur=1;           %maturity at which we want to fit the data. If matur=5, the fifth maturity in the file is chosen.
OptAlg="CMA"; %PatternSearch GeneticAlgorithm SimulatedAnnealing MultiStart
SimPoints=false;
M=50000;           %number of paths to be simulated
iterations=5;     %number of repetitions to be simulated (and then averaged)


%%%%%%%%%%%%%      ORIGINAL DATA MODIFICATIONS     %%%%%%%%%%%%
B(:,2)=B(:,2)/S0;     %normalize strike prices
S0=1;
B(:,1)=B(:,1)/252;    %convert maturities from days to years
times=unique(B(:,1));
T=times(matur);

%B=B(B(:,2)>=0.7 ,:); 
%B(6,:)=[];

B=B(B(:,1)==T,:);     %only keep values of the maturity
%L = T*252*2;          %number of steps in simulations

P=B;
for i=1:size(B,1)
    P(i,3)=european_bs(S0,B(i,2),r,B(i,3),B(i,1),"call");
end


%%%%%%%%%%%%%%%%%%%%%    CALIBRATION      %%%%%%%%%%%%%%%%%%%%%%
x0=0.2;   %parameter starting values [kappa,nubar,nu0,rho,chi]

%%{
sigma=Optimizer(B,x0);
%}


%PlotterPrice(sigma,S0,r,T,P)

PlotterVol(sigma,B)

R=zeros(size(B,1),1);
for i=1:size(B,1)
    R(i)=european_bs(S0,B(i,2),r,sigma,B(i,1),"call");
end

beep
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



%%%%%%%%%%%%%      DEFINE MINIMIZATION PROCEDURE      %%%%%%%%%%%%%%%
function optimvars=Optimizer(B,x0)


fun=@(var)Constcal(var(1),B); %function to be optimized.
%Hestoncal(kappa,nubar,nu0,rho,chi,S0,P,r)

    optimvars=purecmaes(fun,x0);
    

%print optimization output
fprintf('sigma=%f\nerror=%f\n',[optimvars(1),Constcal(optimvars(1),B)])
end



%%%  CALCULATE LEAST SQUARES ERROR BETWEEN MODEL AND DATA  IMPLIED VOL  %%%
function Total_Error=Constcal(sigma,B)

LS=0;
for i=1:size(B,1)
    %LS is the least squares error
    %Function sigmaSABR outputs the error between model and data implied
    %volatilities using the original Hagan formula
    %LS=LS+((P(i,3)-HestonPrice(P(i,2),P(i,1),S0,r,kappa,nubar,nu0,rho,chi,"price"))^2);
    LS=LS+(B(i,3)-sigma)^2*(1-abs(B(i,2)-1))^2;
    %sigmaSABR(alpha,rho,nu,beta,K,f,T)
end
Total_Error=LS;

if sigma<0
   Total_Error=1000; 
end
end




%%%%%%%%%%%%%%%%%%%%%    PLOT OPTIMIZATION RESULTS    %%%%%%%%%%%%%%%%%%%%
function PlotterPrice(sigma,S0,r,T,P)
figure

scatter(P(:,2),P(:,3),100,'x','LineWidth',1.5);


hold on;
 ConstVol=@(K)european_bs(S0,K,r,sigma,T,"call");

 fplot(ConstVol,[0.4,1.6],'LineWidth',1.5)
     xlim([0.4,1.6])
set(gca,'fontsize',12)
xlabel('K/S_0');
ylabel('Option Price (€)')
grid on;
    box on;
pbaspect([1.5 1 1])
end



function PlotterVol(sigma,B)
figure

scatter(B(:,2),B(:,3),100,'x','LineWidth',1.5);
hold on;

 ConstVol=@(K)K*0+sigma;

 fplot(ConstVol,[0.4,1.6],'LineWidth',1.5)
     xlim([0.4,1.6])
     
     
grid on;
set(gca,'fontsize',12)
xlabel('K/S_0');
ylabel('\sigma_{imp}(yr^{-1/2})')
grid on;
box on;
pbaspect([1.5 1 1])
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