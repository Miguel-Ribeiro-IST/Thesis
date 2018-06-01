%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   FOR COMMENTS ON THE CODE, CHECK SIMILAR FILE SABR_Static_Fit.m   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
A = importdata('Data_BNPP_2.txt','\t',1);
B=A.data(:,:);


S0=17099.4;
r = 0;
matur=4;           %maturity until which we want to fit the data.
iterations2=1;
%If matur=5, all maturities until the fifth maturity in the file are chosen.


%beta=0.5;
OptAlg="CMA"; %PatternSearch GeneticAlgorithm SimulatedAnnealing MultiStart
SimPoints=false;
M=100000;
iterations=3;


B(:,2)=B(:,2)/S0;
S0=1;
B(:,1)=B(:,1)/252;
times=unique(B(:,1));
B=B(B(:,1)<=times(matur),:);

x0=[0.2,-0.2,0.8,1,1,0.75];
lb = [0,-1,0,0,0,0];
ub = [5,1,5,100,100,1];

optimvars=IRLSOptim(S0,B,r,x0,OptAlg,lb,ub,iterations2);
alpha=optimvars(1);
rho0=optimvars(2);
nu0=optimvars(3);
a=optimvars(4);
b=optimvars(5);
beta=optimvars(6);


Plotter(alpha,rho0,nu0,a,b,beta,S0,r,B,M,iterations,matur,SimPoints,OptAlg)
beep

function optimvars=IRLSOptim(S0,B,r,x0,OptAlg,lb,ub,iter)
weights=1/size(B,1)*ones(size(B,1),1);
for i=1:iter
    vars=Optimizer(S0,B,r,x0,OptAlg,lb,ub,weights);
    weights=Reweight(vars(1),vars(2),vars(3),vars(4),vars(5),vars(6),S0,B,r);
end
optimvars=vars;
f=SABRcal(optimvars(1),optimvars(2),optimvars(3),optimvars(4),optimvars(5),optimvars(6),S0,B,r,ones(size(B,1),1));
%print optimization output
fprintf(strcat(strcat(strcat("alpha=",strcat(num2str(optimvars(1)),strcat(",    rho0=",strcat(num2str(optimvars(2)),strcat(",    nu0=",strcat(num2str(optimvars(3)),strcat(",    a=",strcat(num2str(optimvars(4)),strcat(",    b=",strcat(num2str(optimvars(5)),strcat(",    beta=",strcat(num2str(optimvars(6)),strcat("\nerror=",strcat(num2str(f),",    ")))))))))))))),OptAlg),"\n"));

end


function weights=Reweight(alpha,rho,nu,a,b,beta,S0,B,r)
errors=ones(size(B,1),1);
for i=1:size(B,1)
    errors(i)=(B(i,3)-sigmaSABR(alpha,rho,nu,a,b,beta,B(i,2),S0*exp(r*B(i,1)),B(i,1)))^2;
end

m=median(errors);
for i=1:size(B,1)
    if errors(i)>6*m
    weights(i)=0;
    else
        weights(i)=(1-(errors(i)/(6*m))^2)^2;
    end
end
weights=weights/sum(weights);

end


function optimvars=Optimizer(S0,B,r,x0,OptAlg,lb,ub,weights)
fun=@(var)SABRcal(var(1),var(2),var(3),var(4),var(5),var(6),S0,B,r,weights);

if OptAlg=="MultiStart"
    rng default % For reproducibility
    opts = optimoptions(@fmincon,'Display','off','Algorithm','sqp');
    %problem = createOptimProblem('fmincon','objective',fun,'x0',x0,'lb',lb,'ub',ub,'options',opts,'nonlcon',@Feller_Condition);
    
    problem = createOptimProblem('fmincon','objective',fun,'x0',x0,'lb',lb,'ub',ub,'options',opts);
    
    ms = MultiStart('UseParallel',true,'StartPointsToRun','bounds-ineqs','Display','off');
    optimvars = run(ms,problem,2);
    %ms = GlobalSearch('StartPointsToRun','bounds-ineqs','Display','off');
    %[optimvars,f] = run(ms,problem);

    elseif OptAlg=="CMA"
    optimvars=purecmaes(fun,x0);
    
end

end


function Plotter(alpha,rho0,nu0,a,b,beta,S0,r,B,M,iterations,matur,SimPoints,OptAlg)
figure
times=unique(B(:,1));
for iter=1:matur
    ax(iter) = subplot(2,ceil(matur/2),iter);
    T=times(iter);
    
    C=B(B(:,1)==T,2:3);
    
    SABRVol=@(K)sigmaSABR(alpha,rho0,nu0,a,b,beta,K,S0*exp(r*T),T);
    
    if SimPoints
    for i=1:size(C,1)
        Volatility(i)= Pricer(alpha,rho0,nu0,a,b,beta,S0,r,T,M,T*252*2,C(i,1),iterations,"vol");
    end
    scatter(ax(iter),C(:,1),Volatility(:),'x');
    hold on;
    end
    scatter(ax(iter),C(:,1),C(:,2),'.');
    hold on;

    fplot(ax(iter),SABRVol,[min(C(:,1)) max(C(:,1))])
    xlim([0.4,1.6])
    ylim([0.2,1])
    hold on;
    title(ax(iter),strcat(strcat(strcat(num2str(T*252)," days  ("),num2str(T*252/21))," months)"))
    clear Volatility

end

text1=strcat(strcat(strcat(num2str(times(1)*252)," days  ("),num2str(times(1)*252/21))," months)");
if SimPoints
    vars1=strcat(strcat(strcat(strcat(strcat("\beta=",num2str(beta)),strcat(",  paths=",num2str(M))),strcat(",  iterations=",num2str(iterations))),",  method="),OptAlg);
else
    vars1=strcat(strcat("\beta=",num2str(beta)),strcat(",  method=",OptAlg));
end
    text2=strcat(strcat(strcat(num2str(times(2)*252)," days  ("),num2str(times(2)*252/21))," months)");
vars2=strcat(strcat(strcat(strcat(strcat(strcat(strcat("\alpha=",num2str(alpha)),strcat(",  \rho_0=",num2str(rho0))),strcat(",  \nu_0=",num2str(nu0))),",  a="),num2str(a)),",  b="),num2str(b));
title(ax(1),{vars1,text1})
title(ax(2),{vars2,text2})
end


function Total_Error=SABRcal(alpha,rho,nu,a,b,beta,S0,B,r,weights)

if alpha<0 || alpha>5 || rho<-1 || rho>1 || nu<0 || nu>5 || a<0 || a>100 || b<0 || b>100 || beta<0 || beta>1
    Total_Error=10000;
else
    
LS=0;
for i=1:size(B,1)
    LS=LS+weights(i).*((B(i,3)-sigmaSABR(alpha,rho,nu,a,b,beta,B(i,2),S0*exp(r*B(i,1)),B(i,1)))/B(i,3))^2;
end
Total_Error=LS;
end
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

parfor iter=1:iterations
    S = S0*ones(M,1);
    sigma=alpha*ones(M,1);
    
    for k = 1:L
        rho=rho0*exp(-a*dt*(k-1));
        nu=nu0*exp(-b*dt*(k-1));
        Z1=randn(M,1);
        Z2=rho*Z1+sqrt(1-rho^2)*randn(M,1);
        S(:)=S(:).*(1+r*dt)+exp(-r*(T-dt*k)*(1-beta)).*sigma(:).*S(:).^beta.*sqrt(dt).*Z1+beta/2*exp(-2*r*(T-dt*k)*(1-beta))*sigma(:).^2.*S(:).^(2*beta-1)*dt.*(Z1.^2-1)
        sigma(:)=sigma(:).*(1+nu*sqrt(dt).*Z2+nu^2/2*dt*(Z2.^2-1));
    end
    
    Y=zeros(M,1);
    for i=1:M
        Y(i) = max(S(i)-K,0);
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