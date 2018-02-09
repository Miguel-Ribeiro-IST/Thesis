clear;
S0=177.12;
r=0.06;

A = importdata('Options_AAPL.txt',';',1);
B=A.data(:,:);
Time=519;
sigma=0.2;
B=B(B(:,1)==Time,2:3);
B=B(abs(B(:,1))<sqrt(Time/252)*sigma*S0,:);
B(:,2)=smooth(B(:,1),B(:,2))

for i=1:size(B,1)
euro=@(sigma)european_bs(S0,B(i,1)+S0,r,sigma,Time/252,'put')-B(i,2);
C(i)=fzero(euro,0.5);
end
%scatter3(B(:,1),B(:,2),C);
%vol=fit([B(:,1),B(:,2)],C','poly22');
%vol(100,100)
%plot(vol,[B(:,1),B(:,2)],C')

F = @(var,x)var(1)+var(2)*x+var(3)*x.*x;
data0 = [0 0 0];
opts=optimset('Display','off');
[var] = lsqcurvefit(F,data0,B(:,1),C',[],[],opts);
vol=@(x)var(1)+var(2)*x+var(3)*x.*x;
scatter(B(:,1),C,'filled');
hold on;
fplot(vol);


%fsurf(@(x,y)2*vol(x,y),[0 600 -200 200]);
%fsurf(vol,[-200 200 0 600])
%plot(F(var,x,y),[B(:,1),B(:,2)],C')




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