clear;
A = importdata('Data_BNPP.txt','\t',1);
B1=A.data(:,:);
S01=17099.4;
r1 = 0.06;
B1(:,2)=B1(:,2)/S01;
S01=1;
B1(:,1)=B1(:,1)/252;
P1=B1;
for i=1:size(B1,1)
    P1(i,3)=european_bs(S01,B1(i,2),r1,B1(i,3),B1(i,1),"call");
end
times1=unique(B1(:,1));


matur1=4;
MinT=min(B1(:,1));
MaxT=max(B1(:,1));
MinK=0.5;
MaxK=1.70;
dT=10.5/252;
dK=0.05*S01;
PriceVol="vol";
iterations=10;
M1=10000;
sigmamax=5;

interpol=Dupire(S01,r1,B1,MinT,MaxT,dT,MinK,MaxK,dK);

tic
figure
for iter=1:matur1
    ax(iter) = subplot(2,ceil(matur1/2),iter);
    T1 = times1(iter);
    L1 = T1*252*2;
    
    if PriceVol=="price"
        C=P1(P1(:,1)==T1,2:3);
    else
        C=B1(B1(:,1)==T1,2:3);
    end
    
    for i=1:size(C,1)
        Euro(i)=Pricer(S01,r1,T1,M1,L1,C(i,1),iterations,PriceVol,sigmamax,interpol);
    end
    
    scatter(ax(iter),C(:,1),C(:,2),'.');
    hold on;
    scatter(ax(iter),C(:,1),Euro(:),'x');
    hold on;
    title(ax(iter),times1(iter)*252)
    clear Euro
end

Timer(0,toc)
beep



function Euro_final=Pricer(S0,r,T,M,L,K,iterations,PriceVol,sigmamax,interpol)
dt = T/L;

parfor iter=1:iterations
    S = S0*ones(M,1);
    sigma=interpol(0,S0)*ones(M,1);
    
    for k = 1:L
        S(:)=S(:)+S(:)*r*dt+sqrt(dt)*sigma(:).*S(:).*randn(M,1);
        
        for i=1:M
           sigma(i)=min(interpol(k*dt,S(i)),sigmamax);
        end
    end
    Y=zeros(M,1);
    for i=1:M
        Y(i) = max(S(i)-K,0);
    end
    
    if PriceVol=="price"
        Euro(iter)=exp(-r*T)*mean(Y(:));
    else
        euro=@(sigma)european_bs(S0,K,r,sigma,T,"call")-exp(-r*T)*mean(Y(:));
        Euro(iter)=fzero(euro,0.25);
    end
end

Euro_final=mean(Euro);
end

function interp=Dupire(S0,r,B,MinT,MaxT,dT,MinK,MaxK,dK)
B=B(B(:,2)<=MaxK & B(:,2)>=MinK & B(:,1)>=MinT & B(:,1)<=MaxT,:);
Time=MinT:2*dT:MaxT;
Strike=MinK:2*dK:MaxK;
[X,Y] = meshgrid(Time,Strike);
S=scatteredInterpolant(B(:,1),B(:,2),B(:,3),'linear','none');
SgradK=(S(X,Y+dK)-S(X,Y-dK))/(2*dK);
Sgrad2K=(S(X,Y+dK)+S(X,Y-dK)-2*S(X,Y))/(dK^2);
SgradT=(S(X+dT,Y)-S(X,Y))/(dT);
SXY=S(X,Y);

d1=(log(S0./Y)+(r+0.5*SXY.^2).*X)./(SXY.*sqrt(X));
vol=sqrt((SXY.^2+2*X.*SXY.*SgradT+2*r*Y.*X.*SXY.*SgradK)./((1+Y.*d1.*sqrt(X).*SgradK).^2+Y.^2.*X.*SXY.*(Sgrad2K-d1.*(SgradK).^2.*sqrt(X))));


V=[];
for i=1:size(vol,1)
   for j=1:size(vol,2)
       if ~isnan(vol(i,j))
           V=[V;[X(i,j),Y(i,j),vol(i,j)]];
       end
   end
end
interp=scatteredInterpolant(V(:,1),V(:,2),V(:,3),'linear','linear');
end

function euro=european_bs(S0,K,r,sigma0,T,putcall)
d1 = (log(S0/K) + (r + 0.5*sigma0^2)*T)/(sigma0*sqrt(T));
d2 = d1 - sigma0*sqrt(T);
N1 = normcdf(d1);
N2 = normcdf(d2);
if putcall=="call"
    euro = S0*N1 - K*exp(-r*T)*N2;
elseif putcall=="put"
    euro = S0*N1 - K*exp(-r*T)*N2 + K*exp(-r*T) - S0;
end
end


function Timer(error,time)
if time>3600
    time=strcat(strcat(strcat(num2str(floor(time/3600)),"hrs,"),num2str(floor((time-floor(time/3600)*3600)/60))),"min");
elseif time>60
    time=strcat(strcat(strcat(num2str(floor(time/60)),"min,"),num2str(floor(time-floor(time/60)*60))),"sec");
else
    time=strcat(num2str(floor(time)),"sec");
end
format shortg;
c = clock;
if error==0
disp(strcat(time,strcat("   ",strcat(num2str(c(4),'%02.f'),strcat(":",num2str(c(5),'%02.f'))))))
fprintf('____________________________________________________\n\n')
else
disp(strcat("error=",num2str(error),strcat("   ",strcat(strcat(time,strcat("   ",strcat(num2str(c(4),'%02.f'),strcat(":",num2str(c(5),'%02.f')))))))))
fprintf('____________________________________________________\n\n')
end
end