clear; %clear all variables. This prevents multiple simulations from interfering with one another
%Import data from file "Data_BNPP.txt" and assign it to matrix B. The file must be in the same folder
%1st column - maturities, 2nd column - strikes, 3rd column - implied volatilities
A = importdata('Data_BNPP.txt','\t',1);
B=A.data(:,:);


%%%%%%%%%%%%%%%%%%%%  INPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%
tic
S0=17099.4;        %initial stock price
r = 0;          %risk-free rate. Forward prices in data file assumed r=0.06
matur=4;           %maturity until which we want to fit the data.
sigmamax=1.5;        %maximum value the local volatility can take
M=10000;           %number of paths to be simulated
aver=10;
%L=T*252*2


%%%%%%%%%%%%%      ORIGINAL DATA MODIFICATIONS     %%%%%%%%%%%%
B(:,2)=B(:,2)/S0;     %normalize strike prices
S0=1;
B(:,1)=B(:,1)/252;    %convert maturities from days to years
times=unique(B(:,1));
B=B(B(:,1)<=times(matur),:);  %restrict data to selected maturity



%%%%%%%%%%   IMPLIED VOLATILITY INTERPOLATION MESH PARAMETERS   %%%%%%%%%
MinT=min(B(:,1));   %minimum value for maturity in the mesh grid
MaxT=max(B(:,1));   %maximum value for maturity in the mesh grid
MinK=0.4;           %minimum value for strike in the mesh grid
MaxK=1.6;          %minimum value for strike in the mesh grid
dT=10.5/252;        %mesh grid size w.r.t. maturity
dK=0.05*S0;         %mesh grid size w.r.t. strike


interpol=Dupire(S0,r,B,MinT,MaxT,dT,MinK,MaxK,dK);

%Plotter(S0,r,B,M,matur,sigmamax,interpol,aver)
tab=Plotter3D(interpol,sigmamax,S0,r,B,M,aver,matur);
%tab=Printer(sigmamax,interpol,M,B,S0,r,aver);
openvar('tab')

beep
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%       FUNCTIONS       %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%S0=starting stock price
%r=risk-free rate
%L=number of simulation time steps
%M=number of simulated paths
%B=input data file
%MinT=minimum value for maturity in the mesh grid
%MaxT=maximum value for maturity in the mesh grid
%MinK=minimum value for strike in the mesh grid
%MaxK=minimum value for strike in the mesh grid
%dT=mesh grid size w.r.t. maturity
%dK=mesh grid size w.r.t. strike
%interpol=MATLAB object corresponding to the interpolated surface of implied vols
%sigmamax=maximum value the local volatility can take

%%%%%%%%%%%%%%     OBTAIN THE LOCAL VOLATILITY SURFACE    %%%%%%%%%%
%%Find function in page 844 of file github.com/Miguel-Ribeiro-IST/Thesis/blob/master/References/Wilmott.pdf
function interp=Dupire(S0,r,B,MinT,MaxT,dT,MinK,MaxK,dK)
B=B(B(:,2)<=MaxK & B(:,2)>=MinK & B(:,1)>=MinT & B(:,1)<=MaxT,:);  %Restrict implied volatility data to the selected interpolation limits

%Produce mesh w.r.t. time and strike
Time=MinT:2*dT:MaxT;
Strike=MinK:2*dK:MaxK;
[X,Y] = meshgrid(Time,Strike);

%Generate interpolating surface of the implied volatilities (linear interpolation, no extrapolation)
S=scatteredInterpolant(B(:,1),B(:,2),B(:,3),'linear','linear');

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
interp=scatteredInterpolant(V(:,1),V(:,2),V(:,3),'linear','linear');
end



%%%%%%%%%%%%%%%%%%%%%    PLOT OPTIMIZATION RESULTS    %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%    PLOT MATURITIES IN DIFFERENT FIGURES  %%%%%%%%%%%%
function Plotter(S0,r,B,M,matur,sigmamax,interpol,aver)
times=unique(B(:,1));   %array with all maturity dates

for iter=1:matur
    figure
    
    T = times(iter);
    C=B(B(:,1)==T,2:3);  %Implied volatilities for the selected maturity
        
    %Plot original data points
    scatter(C(:,1),C(:,2),100,'x','LineWidth',1.5);
    hold on;
    
    %Calculate the implied volatility for an array of Ks and connect the dots to generate a function
    K=0.4:0.01:1.6;
    V=zeros(aver,size(K,2));
    parfor i=1:aver
    V(i,:)=Pricer(S0,r,T,M,T*252*2,K',"vol",sigmamax,interpol);
    end
    V=mean(V);
    plot(K,V,'-','LineWidth',1.5)
    
    %Plot options
    xlim([0.5,1.6])
    ylim([0,1])
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



%%%%%%%%%%%%%%%%%%%%%    PLOT OPTIMIZATION RESULTS    %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%    PLOT MULTIPLE MATURITIES IN A SINGLE FIGURE  %%%%%%%%%%%%
function MultiPlotter(S0,r,B,M,matur,sigmamax,interpol)
figure
times=unique(B(:,1));
for iter=1:matur
    ax(iter) = subplot(2,ceil(matur/2),iter);
    T = times(iter);
    L = T*252*2;
    
    C=B(B(:,1)==T,2:3);
    
    X=0.4:0.01:1.6;
    Y=Pricer(S0,r,T,M,L,X',"vol",sigmamax,interpol);
    
    %   Volatility=Pricer(S0,r,T,M,L,C(:,1),"vol",sigmamax,interpol);
    
    scatter(ax(iter),C(:,1),C(:,2),'.');
    hold on;
    % scatter(ax(iter),C(:,1),Volatility(:),'x');
    plot(ax(iter),X,Y,'-')
    hold on;
    xlim([0.4,1.6])
    ylim([0.2,1])
    title(ax(iter),strcat(strcat(strcat(num2str(T*252)," days  ("),num2str(T*252/21))," months)"))
    clear Volatility
end
text1=strcat(strcat(strcat(num2str(times(1)*252)," days  ("),num2str(times(1)*252/21))," months)");
vars1=strcat(strcat("sigmamax=",num2str(sigmamax)),strcat(",  paths=",num2str(M)));
title(ax(1),{vars1,text1})
end



%%%%%%%%%%%%%%    PLOT OPTIMIZATION RESULTS IN A SURFACE   %%%%%%%%%%%%%%%
function tab=Plotter3D(interpol,sigmamax,S0,r,B,M,aver,matur)
    figure
     times=unique(B(:,1));   %array with all maturity dates
     
    %Plot original data points
    scatter3(B(:,2),B(:,1),B(:,3),30,'LineWidth',1.75,'MarkerEdgeColor','k','MarkerFaceColor',[0.3010    0.7450    0.9330]);
    hold on;
    
    %set variables to plot the 2D-functions (not the surface)
    K2=0.4:0.01:1.6; %!DO NOT CHANGE!
    DV2=zeros(size(times,1),size(K2,2));
    
    %Plot implied volatility surface
    DupireVol=@(K,T)Pricer(S0,r,T,M,T*252*2,K',"vol",sigmamax,interpol);
    [K,T] = meshgrid(K2,0.5/12:0.5/12:0.5+0.5/12); %create the grid to be evaluated
    DV_tmp=zeros(aver,size(K,2)); %matrix to be averaged (reducing noise) to generate the values
    f=1; %auxiliary variable
    for i=1:size(K,1)
        parfor j=1:aver
            DV_tmp(j,:)=DupireVol(K(i,:),T(i,1));
        end
        mn=mean(DV_tmp,1)';
        if i==2 || i==4 || i==6 || i==12
            DV2(f,:)=mn;
            DV2max(f,:)=quantile(DV_tmp,0.975,1);
            DV2min(f,:)=quantile(DV_tmp,0.025,1);
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
    caxis([0 1])
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
caxis([0 1])
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




for iter=1:matur
    figure
    
    T = times(iter);
    C=B(B(:,1)==T,2:3);  %Implied volatilities for the selected maturity
        
    %Plot original data points
    scatter(C(:,1),C(:,2),100,[0    0.1470    0.6410],'x','LineWidth',1.5);
    hold on;
    
    plot(K2,DV2(iter,:),'-.','LineWidth',2,'Color',[0.0010    0.60    0.8330]);
    
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
    %lgd=legend({'Market Data','Simulated Function'},'Location','northeast','FontSize',11);
    %title(lgd,strcat(strcat("T=",num2str(T*252))," days"))
    
    
    h = get(gca,'Children');
    lgd=legend([h(3) h(2) h(1)],{'Market Data','Simulated Function (mean)','95% Confidence Band'},'Location','northeast','FontSize',11);
    title(lgd,strcat(strcat("T=",num2str(T*252))," days"))
    set(gca,'Children',[h(3) h(2) h(1)])
    
    clear Volatility
end




MKTVols=B(:,3);   %Market implied volatilities
MKTPrices=european_bs(S0,B(:,2),r,B(:,3),B(:,1),"call");  %Market (converted) prices

DV3=DV2(:,[11,36,51,61,71,86,111]);

Vols=[]; Prices=[];
K3=unique(B(:,2));
for i=1:matur
    Vols=[Vols;DV3(i,:)'];
    DP=european_bs(S0,K3,r,DV3(i,:)',times(i),"call");
    Prices=[Prices;DP];
end

%Output table
tab=[B(:,1)*252,B(:,2),MKTVols,Vols,abs(MKTVols-Vols)./MKTVols*100,MKTPrices,Prices,abs(MKTPrices-Prices)./MKTPrices*100];
format short
end



%%% CALCULATE MONTE CARLO PRICE/IMPLIED VOLATILITY OF EUROPEAN OPTION UNDER SABR %%%
function Result=Pricer(S0,r,T,M,L,C,PriceVol,sigmamax,interpol)
dt = T/L;      %time steps
N=size(C,1);

%%NOTE:
%%%%We don't need to simulate an entire matrix of stock prices for all time steps and for all paths.
%%%%We just need to simulate one vector of stock prices for all paths and update it at each time step
S = S0*ones(M,1);                  %define initial vector of stock prices
sigma=interpol(0,S0)*ones(M,1);    %define the initial vector of volatilities

for k = 1:L
    S(:)=S(:)+S(:)*r*dt+sqrt(dt)*sigma(:).*S(:).*randn(M,1);   %GBM formula
    
    for i=1:M
        %At each step, calculate the new local volatility value for each path (maximized by threshold "sigmamax")
        sigma(i)=min(interpol(k*dt,S(i)),sigmamax);
    end
end

Y=zeros(M,N);
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
        
        if ~isnan(res)
            Result(j)=res;
        else
            Result(j)=0;
        end
    end
end

end



%%% PRINT A TABLE WITH MODEL/MARKET IMPLIED VOL/PRICES AND REL. ERRORS %%%
function tab=Printer(sigmamax,interpol,M,B,S0,r,aver)
format longG      %Change format for maximum precision
MKTVols=B(:,3);   %Market implied volatilities
MKTPrices=european_bs(S0,B(:,2),r,B(:,3),B(:,1),"call");  %Market (converted) prices

Vols=[]; Prices=[];
    DupireVol=@(K,T)Pricer(S0,r,T,M,T*252*2,K',"vol",sigmamax,interpol);    
    times=unique(B(:,1));   %array with all maturity dates
    K2=unique(B(:,2));
    for i=1:size(times,1)
        DV_tmp=zeros(aver,size(K2,1));
        for f=1:aver
           DV_tmp(f,:)=DupireVol(K2',times(i));
        end
        DV=mean(DV_tmp)';
        DP=european_bs(S0,K2,r,DV,times(i),"call");
       Vols=[Vols;DV];
       Prices=[Prices;DP];
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