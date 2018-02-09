clear;
S0 = 177.12;
K = 177.12+10.38;
T = 4/252;
r = 0.06;
marketeuro=9.35;
putcall='put';
euro=@(sigma)european_bs(S0,K,r,sigma,T,putcall)-marketeuro;
vol=fzero(euro,0.5)



function euro=european_bs(S0,K,r,sigma0,T,putcall)
d1 = (log(S0/K) + (r + 0.5*sigma0^2)*T)/(sigma0*sqrt(T));
d2 = d1 - sigma0*sqrt(T);
N1 = normcdf(d1);
N2 = normcdf(d2);
if putcall=='put'
    euro = S0*N1 - K*exp(-r*T)*N2 + K*exp(-r*T) - S0;
else
    euro = S0*N1 - K*exp(-r*T)*N2;
end
end