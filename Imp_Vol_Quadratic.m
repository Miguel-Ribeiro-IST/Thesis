clear;
S0=167.65;
r=0.06;

A = importdata('AAPL.txt','\t',1);
B=A.data(:,:);

B=B(B(:,1)==35,:);

%%{
for i=1:size(B,1)
euro=@(sigma)european_bs(S0,B(i,2),r,sigma,B(i,1)./252,'call')-B(i,3);
B(i,3)=fzero(euro,0.5);
end
%}

a=0.24;
b=2.1;
c=190;

F = @(var,x)var(1)+var(2)*((x-var(3))/var(3)).^2;
data0 = [a b c];
opts=optimset('Display','off');
[var] = lsqcurvefit(F,data0,B(:,2),B(:,3),[0.1 0.1 165],[0.35 100 200],opts)
vol=@(x)var(1)+var(2)*((x-var(3))/var(3)).^2;
%a1=0.2322;b1=   17.5750;c1=  178.1727;
var2=[0.2103   19.0574  183.3728]
vol2=@(x)var2(1)+var2(2)*((x-var2(3))/var2(3)).^2;

scatter(B(:,2),B(:,3));
hold on;
fplot(vol, [150 185]);
hold on;
fplot(vol2, [150 185]);


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