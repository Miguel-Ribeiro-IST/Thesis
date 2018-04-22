clear;
S0 = 1;
T = 1;
r = .6;
sigma = 0.2;
L = 100; % number of time intervals
dt = T/L;
M = 1; % number of asset paths

rho=1.;
a=0.2;
b=5;
c=0.2;

Y = zeros(M,L);
S = S0*ones(M,L+1); % asset paths
D = S0*ones(M,L+1);
var=sigma^2*ones(M,L+1);
for k = 2:L+1
    X1=randn(M,1);
    Z=randn(M,1);
    X2=rho*X1+sqrt(1-rho^2)*Z;
S(:,k)=S(:,k-1)+r*dt*S(:,k-1)+sqrt(var(:,k-1)).*S(:,k-1)*sqrt(dt).*X1;
D(:,k)=D(:,k-1)+r*dt*D(:,k-1)+sigma*D(:,k-1)*sqrt(dt).*X1;
var(:,k)=var(:,k-1)+(a-b*var(:,k-1))*dt+c*sqrt(var(:,k-1))*sqrt(dt).*X2;
end

plot(S(:,:)')
hold on;
plot(D(:,:)')

