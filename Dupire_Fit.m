clear; %clear all variables. This prevents multiple simulations from interfering with one another
%Import data from file "Data_BNPP.txt" and assign it to matrix B. The file must be in the same folder
%1st column - maturities, 2nd column - strikes, 3rd column - implied volatilities
A = importdata('Data_BNPP.txt','\t',1);
B=A.data(:,:);

%%%%%%%%%%%%%%%%%%%%  INPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%
S0=17099.4;        %initial stock price
r = 0.06;          %risk-free rate. Forward prices in data file assumed r=0.06
matur=4;           %maturity until which we want to fit the data.
                   %If matur=5, all maturities until the fifth maturity in the file are chosen.
M=50000;           %number of paths to be simulated
iterations=10;     %number of repetitions to be simulated (and then averaged)
sigmamax=5;        %maximum value the local volatility can take


%%%%%%%%%%%%%      ORIGINAL DATA MODIFICATIONS     %%%%%%%%%%%%
B(:,2)=B(:,2)/S0;     %normalize strike prices
S0=1;
B(:,1)=B(:,1)/252;    %convert maturities from days to years


%%%%%%%%%%   IMPLIED VOLATILITY INTERPOLATION MESH PARAMETERS   %%%%%%%%%
MinT=min(B(:,1));   %minimum value for maturity in the mesh grid
MaxT=max(B(:,1));   %maximum value for maturity in the mesh grid
MinK=0.5;           %minimum value for strike in the mesh grid
MaxK=1.70;          %minimum value for strike in the mesh grid
dT=10.5/252;        %mesh grid size w.r.t. maturity
dK=0.05*S0;         %mesh grid size w.r.t. strike

tic
interpol=Dupire(S0,r,B,MinT,MaxT,dT,MinK,MaxK,dK);

Plotter(S0,r,B,M,iterations,matur,sigmamax,interpol)
toc

%%%%%%%%%%%%%%%%%%%%%    PLOT INTERPOLATION RESULTS    %%%%%%%%%%%%%%%%%%%%
function Plotter(S0,r,B,M,iterations,matur,sigmamax,interpol)
figure
times=unique(B(:,1));
for iter=1:matur
    ax(iter) = subplot(2,ceil(matur/2),iter);
    T = times(iter);
    L = T*252*2;
    
        C=B(B(:,1)==T,2:3);
        
    for i=1:size(C,1)
        Volatility(i)=Pricer(S0,r,T,M,L,C(i,1),iterations,"vol",sigmamax,interpol);
    end
    
    scatter(ax(iter),C(:,1),C(:,2),'.');
    hold on;
    scatter(ax(iter),C(:,1),Volatility(:),'x');
    hold on;
    title(ax(iter),strcat(strcat(strcat(num2str(T*252)," days  ("),num2str(T*252/21))," months)"))
    clear Volatility
end
text1=strcat(strcat(strcat(num2str(times(1)*252)," days  ("),num2str(times(1)*252/21))," months)");
vars1=strcat(strcat(strcat("sigmamax=",num2str(sigmamax)),strcat(",  paths=",num2str(M))),strcat(",  iterations=",num2str(iterations)));
title(ax(1),{vars1,text1})
end


%%% CALCULATE MONTE CARLO PRICE/IMPLIED VOLATILITY OF EUROPEAN OPTION UNDER SABR %%%
function Result_Avg=Pricer(S0,r,T,M,L,K,iterations,PriceVol,sigmamax,interpol)
dt = T/L;      %time steps

%%NOTE:
%%%%We don't need to simulate an entire matrix of stock prices for all time steps and for all paths.
%%%%We just need to simulate one vector of stock prices for all paths and update it at each time step
parfor iter=1:iterations               %Perform the "for" cycle in parallel
    S = S0*ones(M,1);                  %define initial vector of stock prices
    sigma=interpol(0,S0)*ones(M,1);    %define the initial vector of volatilities
    
    for k = 1:L
        S(:)=S(:)+S(:)*r*dt+sqrt(dt)*sigma(:).*S(:).*randn(M,1);   %GBM formula
        
        for i=1:M
           %At each step, calculate the new local volatility value for each path (maximized by threshold "sigmamax")
           sigma(i)=min(interpol(k*dt,S(i)),sigmamax);
        end
    end
    
    Y=zeros(M,1);
    for i=1:M
        Y(i) = max(S(i)-K,0);      %Calculate the payoff of all paths (assuming calls)
    end
    
    if PriceVol=="price"
        Result(iter)=exp(-r*T)*mean(Y(:)); %Output the discounted expected payoff
    else
        volatility=@(sigma)european_bs(S0,K,r,sigma,T,"call")-exp(-r*T)*mean(Y(:));
        Result(iter)=fzero(volatility,0.25);   %Calculate the expected implied volatility
    end
end

Result_Avg=mean(Result);
end


%%%%%%%%%%%%%%     OBTAIN THE LOCAL VOLATILITY SURFACE    %%%%%%%%%%
%%Find function in page 844 of file github.com/Miguel-Ribeiro-IST/Thesis/blob/master/References/Wilmott.pdf
function interp=Dupire(S0,r,B,MinT,MaxT,dT,MinK,MaxK,dK)
B=B(B(:,2)<=MaxK & B(:,2)>=MinK & B(:,1)>=MinT & B(:,1)<=MaxT,:);  %Restrict implied volatility data to the selected interpolation limits

%Produce mesh w.r.t. time and strike
Time=MinT:2*dT:MaxT;
Strike=MinK:2*dK:MaxK;
[X,Y] = meshgrid(Time,Strike);

%Generate interpolating surface of the implied volatilities (linear interpolation, no extrapolation)
S=scatteredInterpolant(B(:,1),B(:,2),B(:,3),'linear','none');

%Calculate surface gradients
SgradK=(S(X,Y+dK)-S(X,Y-dK))/(2*dK);             %first-derivative w.r.t. K
Sgrad2K=(S(X,Y+dK)+S(X,Y-dK)-2*S(X,Y))/(dK^2);   %second-derivative w.r.t. K
SgradT=(S(X+dT,Y)-S(X,Y))/(dT);                  %first-derivative w.r.t. T
SXY=S(X,Y);

%Generate local volatility surface
d1=(log(S0./Y)+(r+0.5*SXY.^2).*X)./(SXY.*sqrt(X));
vol=sqrt((SXY.^2+2*X.*SXY.*SgradT+2*r*Y.*X.*SXY.*SgradK)./((1+Y.*d1.*sqrt(X).*SgradK).^2+Y.^2.*X.*SXY.*(Sgrad2K-d1.*(SgradK).^2.*sqrt(X))));

%Generate matrix with the local volatility values on the mesh nodes, to be interpolated
V=[];
for i=1:size(vol,1)
   for j=1:size(vol,2)
       if ~isnan(vol(i,j))                 %some mesh nodes will have no data (vol(i,j)=NaN). These values should be rejected
           V=[V;[X(i,j),Y(i,j),vol(i,j)]];
       end
   end
end

%Generate local volatility surface
interp=scatteredInterpolant(V(:,1),V(:,2),V(:,3),'linear','nearest');
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