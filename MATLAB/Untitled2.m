
matur1=6;
sigma0=0.2;
iterations2=10;
figure

for iter=1:matur1
ax(iter) = subplot(2,matur1/2,iter);
ti=times1(iter);
C=B1(B1(:,1)==ti,2:3);

T1 = times1(iter)/252;
D1=252;
L1 = T1*D1*2;

for i=1:size(C,1)
Euro(i)=Pricer(a1,b1,c1,d1,sigmamax1,S01,r1,T1,D1,M1,L1,C(i,1),iterations2,times1);
end

for i=1:size(C,1)
Euro_Const(i)=Pricer_Const(sigma0,S01,r1,T1,M1,L1,C(i,1),iterations2);
end


scatter(ax(iter),C(:,1),C(:,2));
hold on;
scatter(ax(iter),C(:,1),Euro(:));
hold on;
%scatter(ax(iter),C(:,1),Euro_Const(:));
title(ax(iter),times1(iter)) 
clear Euro Euro_Const
end



function Euro_final=Pricer(a,b,c,d,sigmamax,S0,r,T,D,M,L,K,iterations,times)
dt = T/L;

for iter=1:iterations
S = S0*ones(M,1); % asset paths
sigma=a(1)*ones(M,1);

for k = 2:L+1
S(:)=S(:)+S(:)*r*dt+sqrt(dt)*sigma(:).*S(:).*randn(M,1);

for i=1:size(times,1)
if dt*D*(k-1)<=times(i)
   %sigma=arrayfun(@(x) min(a(i)+b(i)*min(x,0)^2,sigmamax),(S(:,k)-c(i))/c(i));
   sigma=arrayfun(@(x) min(a(i)+(-b(i)*(abs(x-d(i))-(x-d(i)))+c(i)*(abs(x-d(i))+(x-d(i)))).*(x-d(i)),sigmamax),(S(:)));
   %sigma=arrayfun(@(x) a(i)+b(i)*x^2,(S(:,k)-c(i))/c(i));
   break;
end
end

end

for i=1:M
Y(i) = max(S(i)-K,0);
end

   euro=@(sigma)european_bs(S0,K,r,sigma,T,'call')-exp(-r*T)*mean(Y(:));
   Euro_Vol(iter)=fzero(euro,0.2);

end

Euro_final=mean(Euro_Vol);
end


function Euro_final_Const=Pricer_Const(sigma,S0,r,T,M,L,K,iterations)
dt = T/L;

for iter=1:iterations
S = S0*ones(M,L+1); % asset paths

for k = 2:L+1
S(:,k)=S(:,k-1)+S(:,k-1)*r*dt+sqrt(dt)*sigma.*S(:,k-1).*randn(M,1);

end

for i=1:M
Y(i) = max(S(i,L+1)-K,0);
end

Euro_Vol(iter)=exp(-r*dt*L)*mean(Y(:));
end

Euro_final_Const=mean(Euro_Vol);
end


function euro=european_bs(S0,K,r,sigma0,T,putcall)
d1 = (log(S0/K) + (r + 0.5*sigma0^2)*T)/(sigma0*sqrt(T));
d2 = d1 - sigma0*sqrt(T);
N1 = normcdf(d1);
N2 = normcdf(d2);
if putcall=='call'
    euro = S0*N1 - K*exp(-r*T)*N2;
elseif putcall=='put'
    euro = S0*N1 - K*exp(-r*T)*N2 + K*exp(-r*T) - S0;
end
end