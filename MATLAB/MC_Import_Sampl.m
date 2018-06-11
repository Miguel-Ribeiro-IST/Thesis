clear;
S0 = 1;
K = 1.5;
T = 1;
r = 0.05;
sigma = 0.3;
L = 252; % number of time intervals
dt = T/L;
M = 100000; % number of asset paths


r1=log(K/S0);
%r1=1/T*log(K*1.5/S0);

% findfunc=@(rval)0.5*erfc(1./(sqrt(2*(exp(sigma^2*T)-1)))*(1-K./(S0*exp(rval*T))))-0.1;
% r1=fzero(findfunc,0.3);


theta=(r-r1)/sigma;

for j=1:10
    
    
    S = S0*ones(M,1); % asset paths
    S1 = S0*ones(M,1); % asset paths
    W=zeros(M,1);
    for k = 2:L+1
        Z1=randn(M/2,1);
        Z=[Z1;-Z1];
        W=W+sqrt(dt)*Z;
        S(:)=S(:).*exp((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*Z);
        S1(:)=S1(:).*exp((r1-0.5*sigma^2)*dt+sigma*sqrt(dt)*Z);
    end
    RN=exp(-0.5*theta^2*T+theta*W);
    
    Y = zeros(M,1);
    Y1 = zeros(M,1);
    % Find payoff Y at expiry.
    for i=1:M
        Y(i) = max(S(i)-K,0);
        Y1(i) = max(S1(i)-K,0);
    end
    
    
    
    simul_simpl(j)=exp(-r*T)*mean(Y(:));
    simul_simpl1(j)=exp(-r*T)*mean(Y1.*RN);
end

[abs(european_bs(S0,K,r,sigma,T,"call")-mean(simul_simpl1))/abs(european_bs(S0,K,r,sigma,T,"call")-mean(simul_simpl)),std(simul_simpl1)/std(simul_simpl)]

beep

function price=european_bs(S0,K,r,sigma,T,putcall)
d1 = (log(S0./K) + (r + 0.5.*sigma.^2).*T)/(sigma.*sqrt(T));
d2 = d1 - sigma.*sqrt(T);
N1 = normcdf(d1);
N2 = normcdf(d2);
if putcall=="call"
    price = S0.*N1 - K.*exp(-r.*T).*N2;
elseif putcall=="put"
    price = S0.*N1 - K.*exp(-r*T).*N2 + K.*exp(-r.*T) - S0;
end
end