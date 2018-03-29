clear;
S0 = 42;
K = 40;
T = 2;
r = 0.06;
sigma = .3;
L = 50; % number of time intervals
dt = T/L;
M = 100; % number of asset paths

rho=.5;
c=0.3;
b=20;
avg=0.3;
a=avg^2*b;

Y = zeros(M,L);
S = S0*ones(M,L+1); % asset paths
var=sigma^2*ones(M,L+1);
for k = 2:L+1
    X1=randn(M,1);
    Z=randn(M,1);
    X2=rho*X1+sqrt(1-rho^2)*Z;
    S(:,k)=S(:,k-1)+r*dt*S(:,k-1)+sqrt(var(:,k-1)).*S(:,k-1)*sqrt(dt).*X1;
    var(:,k)=var(:,k-1)+(a-b*var(:,k-1))*dt+c*sqrt(var(:,k-1))*sqrt(dt).*X2;
end

figure
ax1 = subplot(2,1,1);
plot(ax1,0:dt:T,S(1:2,:)')
axis(ax1,[0 T 0 3*S0])
title(ax1,'Prices') 

ax2 = subplot(2,1,2);
plot(ax2,0:dt:T,sqrt(var(1:2,:))')
axis(ax2,[0 T 0 2*sigma])
title(ax2,'Volatilities') 


% Find payoff Y at expiry.
for i=1:M
Y(i,L) = max(S(i,L+1)-K,0);
end

d1 = (log(S0/K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
N1 = normcdf(d1);
N2 = normcdf(d2);
European_call_BS = S0*N1 - K*exp(-r*T)*N2;

European_call=exp(-r*dt*L)*mean(Y(:,L));

% Find payoff Y at nodes for each time index.
for k = L+1:-1:3
j = 0;
for i=1:M
if S(i,k-1) > K % in-the-money condition
j = j+1;
S1(j) = S(i,k-1); % in-the-money asset price
Y1(j) = exp(-r*dt)*Y(i,k-1); % discounted cash flow
end
end
p = polyfit(S1,Y1,3);
%F = @(var,x)var(1)+var(2)*x+var(3)*x.*x;
%F = @(var,x)var(1)+var(2)*exp(-x/2)+var(3)*exp(-x/2).*(1-x)+var(4)*exp(-x/2).*(1-2*x+0.5*x.*x);
%data0 = [1 1 1 1];
%opts=optimset('Display','off');
%[var] = lsqcurvefit(F,data0,S1.',Y1.',[],[],opts);

for i = 1:M
if S(i,k-1)-K > polyval(p,S(i,k-1)) % early exercise condition
%if K - S(i,k-1) > F(var,S(i,k-1))
Y(i,k-2) = max(S(i,k-1)-K,0);
else
Y(i,k-2) = exp(-r*dt)*Y(i,k-1);
end
end
end
American_call = exp(-r*dt)*mean(Y(:,1));