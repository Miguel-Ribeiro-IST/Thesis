clear;

T=1;
Smax=200;
N=40000;
M=1000;
S0=40;

K=40;
r=0.06;
sigma=0.2;
q=0;

dt=T/N;
dS=Smax/M;

a = -1/(1+r*dt)*0.5*dt*((r-q)*(0:M)-sigma^2*(0:M).^2);
b = 1/(1+r*dt)*(1 - dt*sigma^2*(0:M).^2);
c = 1/(1+r*dt)*(0.5*dt*(sigma^2*(0:M).^2 + (r-q).*(0:M)));

V=zeros(M+1,N+1);
V(:,N+1) = max(K-(0:dS:Smax)',0);
V(1,:) = K;

for i = N:-1:1
for j=2:M
    V(j,i)=a(j)*V(j-1,i+1)+b(j)*V(j,i+1)+c(j)*V(j+1,i+1);
end
V(:,i)=max(K-[0:dS:Smax]',V(:,i));
end

V(S0/dS+1,1)