%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                    %%%%
%%%%   FOR COMMENTS ON THE CODE, CHECK SIMILAR FILE SABR_Static_Fit.m   %%%%
%%%%                                                                    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except alpha rho0 nu0 beta a b;
A = importdata('Data_BNPP.txt','\t',1);
B=A.data(:,:);


%%%%%%%%%%%%%%%%%%%%  INPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%
S0=17099.4;
r = 0;
matur=4;           %maturity until which we want to fit the data. If matur=5, all maturities until the fifth maturity in the file are chosen.
OptAlg="CMA";      %"CMA" or "MultiStart" optimization algorithms


%%%%%%%%%%%%%%%%%%%   MONTE CARLO SIMULATION %%%%%%%%%%%%%%
SimPoints=false;
M=100000;
repetitions=100;
%L=T*252*2


%%%%%%%%%%%%%      ORIGINAL DATA MODIFICATIONS     %%%%%%%%%%%%
B(:,2)=B(:,2)/S0;
S0=1;
B(:,1)=B(:,1)/252;
times=unique(B(:,1));
B=B(B(:,1)<=times(matur),:);

%%%%%%%%%%%%%%%%%%%%%    CALIBRATION      %%%%%%%%%%%%%%%%%%%%%%
%[alpha, rho0, nu0, a, b, beta]
x0 = [0.2, -0.42, 2.45, 1.14,   2.62,   0.75];
lb = [0,   -1,   0,   0,   0,   0];
ub = [5,   1,    5,   25, 25, 1];
%{
optimvars=Optimizer(S0,B,r,x0,OptAlg,lb,ub);
alpha=optimvars(1);
rho0=optimvars(2);
nu0=optimvars(3);
a=optimvars(4);
b=optimvars(5);
beta=optimvars(6);
%}

tic
%%%%%%%%%%%%%%%%%%%    PLOT RESULTS    %%%%%%%%%%%%%%%%%%%%%
%Plotter(alpha,rho0,nu0,a,b,beta,S0,r,B,M,matur,SimPoints)
Plotter_Sim(alpha,rho0,nu0,a,b,beta,S0,r,B,M,matur,repetitions)
Plotter3D(alpha,rho0,nu0,a,b,beta,S0,r,B)

tab=Printer(alpha,rho0,nu0,a,b,beta,B,S0,r);
%openvar('tab')
toc
beep



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%       FUNCTIONS       %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function optimvars=Optimizer(S0,B,r,x0,OptAlg,lb,ub)
fun=@(var)SABRcal(var(1),var(2),var(3),var(4),var(5),var(6),S0,B,r,lb,ub);

if OptAlg=="MultiStart"
    StartPoints=2;
    opts = optimoptions(@fmincon,'Display','off','Algorithm','active-set');
    problem = createOptimProblem('fmincon','objective',fun,'x0',x0,'lb',lb,'ub',ub,'options',opts);
    ms = MultiStart('UseParallel',true,'StartPointsToRun','bounds-ineqs','Display','off');
    [optimvars,f] = run(ms,problem,StartPoints);
    
elseif OptAlg=="CMA"
    optimvars=purecmaes(fun,x0);
    f=SABRcal(optimvars(1),optimvars(2),optimvars(3),optimvars(4),optimvars(5),optimvars(6),S0,B,r,lb,ub);
end

fprintf('alpha=%f\nrho0=%f, nu0=%f\na=%f,     b=%f\nbeta=%f\n%s, error=%f\n\n',[optimvars(1),optimvars(2),optimvars(3),optimvars(4),optimvars(5),optimvars(6),OptAlg,f])
end



function Total_Error=SABRcal(alpha,rho,nu,a,b,beta,S0,B,r,lb,ub)
if alpha<lb(1) || alpha>ub(1) || rho<lb(2) || rho>ub(2) || nu<lb(3) || nu>ub(3) || a<lb(4) || a>ub(4) || b<lb(5) || b>ub(5) || beta<lb(6) || beta>ub(6)
    Total_Error=10000;
else
    
    LS=0;
    for i=1:size(B,1)
        LS=LS+(B(i,3)-sigmaSABR(alpha,rho,nu,a,b,beta,B(i,2),S0*exp(r*B(i,1)),B(i,1)))^2*(1-abs(B(i,2)-1))^2;
    end
    Total_Error=LS;
end
end



%%%% DYNAMIC SABR CLOSED-FORM SOLUTION FOR IMPLIED VOLATILITY %%%%
% For more information check https://github.com/Miguel-Ribeiro-IST/Thesis/blob/master/References/Fernandez_SABR.pdf
function sigma=sigmaSABR(alpha,rho0,nu0,a,b,beta,K,f,T)
w=alpha.^(-1).*f.^(1-beta);
n1=@(T)2.*nu0*rho0./(T.^2.*(a+b).^2).*(exp(-(a+b).*T)-(1-(a+b).*T));
n22=@(T)3*nu0^2*rho0^2./(T.^4*(a+b)^4).*(exp(-2*(a+b).*T)-8*exp(-(a+b)*T)+(7+2*(a+b)*T.*(-3+(a+b)*T)));
v12=@(T)6*nu0^2./(2*b.*T).^3.*(((2*b*T).^2./2-2*b*T+1)-exp(-2*b*T));
v22=@(T)6*nu0^2./(2*b*T).^3.*(2*(exp(-2*b*T)-1)+2*b*T.*(exp(-2*b*T)+1));
A1=@(T)(beta-1)/2+n1(T).*w/2;
A2=@(T)(1-beta)^2/12+(1-beta-n1(T).*w)/4+(4*v12(T)+3*(n22(T)-3*(n1(T)).^2)).*w.^2/24;
B=@(T)1./w.^2.*((1-beta)^2/24+w.*beta.*n1(T)./4+(2*v22(T)-3*n22(T)).*w.^2./24);
sigma=1./w.*(1+A1(T).*log(K./f)+A2(T).*(log(K./f)).^2+B(T).*T);
end



function Plotter(alpha,rho0,nu0,a,b,beta,S0,r,B,M,matur,SimPoints)
times=unique(B(:,1));

for iter=1:matur
    figure
    
    T=times(iter);
    C=B(B(:,1)==T,2:3);
    
    if SimPoints
        Volatility= Pricer(alpha,rho0,nu0,a,b,beta,S0,r,T,M,T*252*2,C,"vol");
        
        scatter(C(:,1),Volatility(:),100,'+','LineWidth',1.5);
        hold on;
    end
    
    scatter(C(:,1),C(:,2),100,'x','LineWidth',1.5);
    hold on;
    
    SABRVol=@(K)sigmaSABR(alpha,rho0,nu0,a,b,beta,K,S0*exp(r*T),T);
    fplot(SABRVol,[0.4,1.6],'LineWidth',1.5)
    hold on;
    
    xlim([0.4,1.6])
    ylim([0,1])
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
    
    clear Volatility

end
end



function Plotter_Sim(alpha,rho0,nu0,a,b,beta,S0,r,B,M,matur,repetitions)
times=unique(B(:,1));   %array with all maturity dates

K2=0.4:0.01:1.6;
Mdl=zeros(matur,size(K2,2));

scatter3(B(:,2),B(:,1),B(:,3),30,'LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.3010    0.7450    0.9330]);
hold on;
SimVol=@(K,T)Pricer(alpha,rho0,nu0,a,b,beta,S0,r,T,M,T*252*2,K',"vol");
[K,T] = meshgrid(K2,0.5/12:0.5/12:0.5+0.5/12); %create the grid to be evaluated
DV_tmp=zeros(repetitions,size(K,2)); %matrix to be averaged (reducing noise) to generate the values
f=1; %auxiliary variable
for i=1:size(K,1)
    parfor j=1:repetitions
        DV_tmp(j,:)=SimVol(K(i,:),T(i,1));
    end
    mn=mean(DV_tmp,1)';
    if i==2 || i==4 || i==6 || i==12
        DV2(f,:)=mn;
        DV2max(f,:)=quantile(DV_tmp,0.9,1);
        DV2min(f,:)=quantile(DV_tmp,0.1,1);
        plot3(K2,ones(1,size(K2,2))*T(i,1),DV2(f,:),'-.','LineWidth',2,'Color',[0.9500    0.200    0.1])
        hold on;
        f=f+1;
    end
    DV(i,:)=mn;
end
s=surf(K,T,DV);
s.EdgeAlpha=0.6;
s.FaceAlpha=0.85;
shading interp
hold on;

%Plot options
xlim([0.4,1.6])
ylim([0.5/12,0.5+0.5/12])
zlim([0,1])
box on;
grid on;
xlabel('K/S_0');
ylabel('T (days)');
zlabel('\sigma_{imp} (yr^{-1/2})')
yticks([1/12,2/12,3/12,4/12,5/12,6/12])
yticklabels({'21','42','63','84','105','126'})
pbaspect([1 1.5 1])
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
contourf(K,T,DV,25)
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
hold on;



for iter=1:matur
    figure
    
    T=times(iter);
    C=B(B(:,1)==T,2:3);  %Implied volatilities for the selected maturity
    
    %Plot original data points
    scatter(C(:,1),C(:,2),100,[0    0.1470    0.6410],'x','LineWidth',1.5);
    hold on;
    
    %Plot implied volatility function under the Heston model
    SABRVol=@(K)sigmaSABR(alpha,rho0,nu0,a,b,beta,K,S0*exp(r*T),T);
    K=(0.4:0.01:1.6);
    for i=1:size(K,2)
        SV(i)=SABRVol(K(i));
    end
    plot(K,SV,'LineWidth',2,'Color',[0.9500    0.2250    0.0580]);
    hold on;
    
    plot(K,DV2(iter,:),'-.','LineWidth',2,'Color',[0.0010    0.60    0.8330]);
    
    hold on;
    K3 = [K2, fliplr(K2)]; % Use ; instead of ,
    inBetween = [DV2max(iter,:), fliplr(DV2min(iter,:))]; % Use ; instead of ,
    fill(K3, inBetween,[0    0.150    0.830],'FaceAlpha',0.2,'EdgeAlpha',0);
    
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
    lgd=legend([h(4) h(3) h(2) h(1)],{'Market Data','Theoretical Function','Simulated Function (mean)','90% Confidence Interval'},'Location','northeast','FontSize',11);
    title(lgd,strcat(strcat("T=",num2str(T*252))," days"))
    set(gca,'Children',[h(4) h(2) h(3) h(1)])
    
    clear Volatility
end

end



function MultiPlotter(alpha,rho0,nu0,a,b,beta,S0,r,B,M,matur,SimPoints,OptAlg)
figure
times=unique(B(:,1));

for iter=1:matur
    
    if matur>1
        ax(iter) = subplot(2,ceil(matur/2),iter);
    else
        ax(iter) = subplot(1,1,1);
    end
    
    T=times(iter);
    C=B(B(:,1)==T,2:3);
    
    
    if SimPoints
        Volatility= Pricer(alpha,rho0,nu0,a,b,beta,S0,r,T,M,T*252*2,C,"vol");
        
        scatter(ax(iter),C(:,1),Volatility(:),'x');
        hold on;
    end
    scatter(ax(iter),C(:,1),C(:,2),'.');
    hold on;
    
    SABRVol=@(K)sigmaSABR(alpha,rho0,nu0,a,b,beta,K,S0*exp(r*T),T);
    fplot(ax(iter),SABRVol,[0.4,1.6])
    hold on;
    title(ax(iter),strcat(strcat(strcat(num2str(T*252)," days  ("),num2str(T*252/21))," months)"))
    clear Volatility
    xlim([0.4,1.6])
    ylim([0.2,1])
end
text1=strcat(strcat(strcat(num2str(times(1)*252)," days  ("),num2str(times(1)*252/21))," months)");
if SimPoints
    vars1=strcat(strcat(strcat(strcat("\beta=",num2str(beta)),strcat(",  paths=",num2str(M))),",  method="),OptAlg);
else
    vars1=strcat(strcat("\beta=",num2str(beta)),strcat(",  method=",OptAlg));
end
text2=strcat(strcat(strcat(num2str(times(2)*252)," days  ("),num2str(times(2)*252/21))," months)");
vars2=strcat(strcat(strcat(strcat(strcat(strcat(strcat("\alpha=",num2str(alpha)),strcat(",  \rho_0=",num2str(rho0))),strcat(",  \nu_0=",num2str(nu0))),",  a="),num2str(a)),",  b="),num2str(b));
title(ax(1),{vars1,text1})
title(ax(2),{vars2,text2})
end



%%%%%%%%%%%%%%    PLOT OPTIMIZATION RESULTS IN A SURFACE   %%%%%%%%%%%%%%%
function Plotter3D(alpha,rho0,nu0,a,b,beta,S0,r,B)
    figure
    
    %Plot original data points
    scatter3(B(:,2),B(:,1),B(:,3),30,'LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.3010    0.7450    0.9330]);
    hold on;
    
    %Plot implied volatility function under the Heston model
    SABRVol=@(K,T)sigmaSABR(alpha,rho0,nu0,a,b,beta,K,S0*exp(r*T),T);
    [K,T] = meshgrid(0.4:0.01:1.6,0.5/12:0.1/12:0.5+0.5/12);
    for i=1:size(K,1)
        for j=1:size(K,2)
            SV(i,j)=SABRVol(K(i,j),T(i,j));
        end
    end
    s=surf(K,T,SV);
    %s.EdgeColor = 'interp';
    s.EdgeAlpha=0.6;
    s.FaceAlpha=0.85;
    shading interp
    hold on;
        
    times=unique(B(:,1));   %array with all maturity dates
    K2=0.4:0.01:1.6;
    for i=1:size(times,1)
        for j=1:size(K2,2)
           SV2(j)=SABRVol(K2(j),times(i));
        end
        plot3(K2,ones(1,size(K2,2))*times(i),SV2,'LineWidth',2,'Color',[0.9500    0.200    0.1])
    end
    
   xlim([0.4,1.6])
    ylim([0.5/12,0.5+0.5/12])
    zlim([0,1])
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
contourf(K,T,SV,25)
pbaspect([1.5 1 1])
xlabel('K/S_0');
ylabel('T (days)');
yticks([1/12,2/12,3/12,4/12,5/12,6/12])
yticklabels({'21','42','63','84','105','126'})
box on;
colorbar;
set(gca,'fontsize',12)
end



function Result=Pricer(alpha,rho0,nu0,a,b,beta,S0,r,T,M,L,C,PriceVol)
dt = T/L;
N=size(C,1);
S = S0*ones(M,1);
sigma=alpha*ones(M,1);

for k = 1:L
    rho=rho0*exp(-a*dt*(k-1));
    nu=nu0*exp(-b*dt*(k-1));
    Z1=randn(M,1);
    Z2=rho*Z1+sqrt(1-rho^2)*randn(M,1);
    S(:)=max(S(:),0).*(1+r*dt)+exp(-r*(T-dt*k)*(1-beta)).*sigma(:).*max(S(:),0).^beta.*sqrt(dt).*Z1+beta/2*exp(-2*r*(T-dt*k)*(1-beta))*sigma(:).^2.*max(S(:),0).^(2*beta-1)*dt.*(Z1.^2-1);
    sigma(:)=sigma(:).*(1+nu*sqrt(dt).*Z2+nu^2/2*dt*(Z2.^2-1));
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



function tab=Printer(alpha,rho0,nu0,a,b,beta,B,S0,r)
format longG
MKTVols=B(:,3);
MKTPrices=european_bs(S0,B(:,2),r,B(:,3),B(:,1),"call");

SABRVol=@(K)sigmaSABR(alpha,rho0,nu0,a,b,beta,K,S0.*exp(r.*B(:,1)),B(:,1));
Vols=SABRVol(B(:,2));
Prices=european_bs(S0,B(:,2),r,Vols,B(:,1),"call");


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