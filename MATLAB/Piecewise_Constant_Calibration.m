clear;
A = importdata('Data_BNPP.txt','\t',1);
B1=A.data(:,:);
S01=17099.4;
r1 = 0.06;


matur1=2;
M1 = 10;
iterations1=3;
PriceVol="price";
sigmamax1=5;
iterations2=iterations1;



B1(:,2)=B1(:,2)/S01;
S01=1;

times1=unique(B1(:,1));
a1=zeros(matur1,1);
b1=zeros(matur1,1);
c1=zeros(matur1,1);
d1=zeros(matur1,1);

P1=B1;
for i=1:size(B1,1)
    P1(i,3)=european_bs(S01,B1(i,2),r1,B1(i,3),B1(i,1)./252,"call");
end

for iter=1:matur1
    tic
    ti1=times1(iter);
    T1 = ti1/252;
    D1=252;
    L1 = T1*D1*2;
    if PriceVol=="price"
        B1tmp=P1(P1(:,1)==times1(iter),2:3);
    else
        B1tmp=B1(B1(:,1)==times1(iter),2:3);
    end
    
    fun = @(var)Localvol(var(1),var(2),var(3),var(4),iter,S01,r1,T1,D1,M1,L1,B1tmp,a1,b1,c1,d1,sigmamax1,iterations1,times1,PriceVol);
    
    
    A = [];b = [];Aeq = [];beq = [];nonlcon=[];
    lb = [0.15,0.05,0,0.9];
    ub = [0.25,0.6,1,1.3];
    x0=[0.24-0.002*ti1, 0.2+0.3/ti1, 0.5/iter, 1.15];
    
    options = optimoptions('patternsearch','Display','off','MaxIter',10000,'UseParallel',true);
    vars=patternsearch(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    disp(vars);
    a1(iter)=vars(1);
    b1(iter)=vars(2);
    c1(iter)=vars(3);
    d1(iter)=vars(4);
    

    Timer(Localvol(var(1),var(2),var(3),var(4),iter,S01,r1,T1,D1,M1,L1,B1tmp,a1,b1,c1,d1,sigmamax1,iterations1,times1,PriceVol),toc);
end




figure
for iter=1:matur1
    ax(iter) = subplot(2,ceil(matur1/2),iter);
    ti=times1(iter);
    T1 = times1(iter)/252;
    D1=252;
    L1 = T1*D1*2;
    
    if PriceVol=="price"
        C=P1(P1(:,1)==ti,2:3);
    else
        C=B1(B1(:,1)==ti,2:3);
    end
    
    for i=1:size(C,1)
        Euro(i)=Pricer(a1,b1,c1,d1,sigmamax1,S01,r1,T1,D1,M1,L1,C(i,1),iterations2,times1,PriceVol);
    end
    
    scatter(ax(iter),C(:,1),C(:,2),'.');
    hold on;
    scatter(ax(iter),C(:,1),Euro(:),'x');
    hold on;
    title(ax(iter),times1(iter))
    clear Euro
end

beep



function Error_final=Localvol(ai,bi,ci,di,matur,S0,r,T,D,M,L,B,a,b,c,d,sigmamax,iterations,times,PriceVol)

dt = T/L;
Error=zeros(iterations,1);

a(matur)=ai;
b(matur)=bi;
c(matur)=ci;
d(matur)=di;

parfor iter=1:iterations
    S = S0*ones(M,1);
    sigma=a(1)*ones(M,1);
    
    for k = 2:L+1
        S(:)=S(:)+S(:)*r*dt+sqrt(dt)*sigma(:).*S(:).*randn(M,1);
        
        for i=1:size(times,1)
            if dt*D*(k-1)<=times(i)
                sigma=arrayfun(@(x) min(a(i)+(-b(i)*(abs(x-d(i))-(x-d(i)))+c(i)*(abs(x-d(i))+(x-d(i)))).*(x-d(i)),sigmamax),(S(:)));
                break;
            end
        end
    end
    
    Y=zeros(M,size(B,1));
    for j=1:size(B,1)
        for i=1:M
            Y(i,j) = max(S(i)-B(j,1),0);
        end
    end
    
    if PriceVol=="price"
        Error(iter)=sum((exp(-r*T)*mean(Y)-B(:,2)').^2);
    else
        Z=exp(-r*T)*mean(Y);
        for i=1:size(B,1)
            euro=@(sigma)european_bs(S0,B(i,1),r,sigma,T,"call")-Z(i);
            Z(i)=fzero(euro,0.25);
        end
        Error(iter)=sum((Z-B(:,2)').^2);
    end
    
end
Error_final=mean(Error);
end



function Euro_final=Pricer(a,b,c,d,sigmamax,S0,r,T,D,M,L,K,iterations,times,PriceVol)
dt = T/L;

for iter=1:iterations
    S = S0*ones(M,1);
    sigma=a(1)*ones(M,1);
    
    for k = 2:L+1
        S(:)=S(:)+S(:)*r*dt+sqrt(dt)*sigma(:).*S(:).*randn(M,1);
        
        for i=1:size(times,1)
            if dt*D*(k-1)<=times(i)
                sigma=arrayfun(@(x) min(a(i)+(-b(i)*(abs(x-d(i))-(x-d(i)))+c(i)*(abs(x-d(i))+(x-d(i)))).*(x-d(i)),sigmamax),(S(:)));
                break;
            end
        end
        
    end
    
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
if time>60
    time=strcat(num2str(floor(time/60)),"min");
else
    time=strcat(num2str(floor(time)),"sec");
end
format shortg;
c = clock;
disp(strcat("error=",num2str(error),strcat("   ",strcat(strcat(time,strcat("   ",strcat(num2str(c(4),'%02.f'),strcat(":",num2str(c(5),'%02.f')))))))))
fprintf('____________________________________________________\n\n')
end