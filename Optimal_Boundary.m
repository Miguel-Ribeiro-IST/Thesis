clear;
tic
S0 = 42;
K = 40;
T = 2;
r = 0.06;
sigma = 0.4;
L = 50; % number of time intervals
dt = T/L;
M = 100000; % number of asset paths

val=0;idx=0;
Y = zeros(M,1);
S = S0*ones(M,L+1); % asset paths
for k = 2:L+1
S(:,k)=S(:,k-1).*exp((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*randn(M,1));
end

%S=[1,1.09,1.08,1.34;1,1.16,1.26,1.54;1,1.22,1.07,1.03;1,0.93,0.97,0.92;1,1.11,1.56,1.52;1,0.76,0.77,0.9;1,0.92,0.84,1.01;1,0.88,1.22,1.34];
%T=3;
%K=1.1;
%L=3;
%M=8;
%dt = T/L;

%plot(S(1:100,:)')

B(L+1)=K;

for i=1:M
Y(i) = max(K - S(i,L+1),0);
end

disc=exp(-r*dt);

for t=L:-1:2
    val=0;
    for b=cat(1,0,S(:,t))'
        if b>=K
            continue
        end
        Y1=[];
        for i=1:M
            if S(i,t)<=b
                Y1(i)=K - S(i,t);
            else
                Y1(i)=Y(i)*disc;
            end
        end
        tp=sum(Y1)/M;
        if tp>val
            val=tp;
            B(t)=b;
        end
        end
        for k=1:M
        if S(k,t)<=B(t)
            Y(k)=K - S(k,t);
        else
            Y(k)=Y(k)*disc;
        end
        end
end


M = 1000000; % number of asset paths

Y = zeros(M,L);
S = S0*ones(M,L+1); % asset paths
for k = 2:L+1
S(:,k)=S(:,k-1).*exp((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*randn(M,1));
end

for t=L:-1:2
        for k=1:M
        if S(k,t)<=B(t)
            Y(k)=K - S(k,t);
        else
            Y(k)=Y(k)*disc;
        end
        end
end
toc
beep
American=val*disc