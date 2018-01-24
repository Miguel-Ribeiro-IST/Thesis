clear;
S0=177.12;
r=0.06;

A = importdata('Options_AAPL.txt',';',1);
B=A.data(:,:);
%{
Time=174;
sigma=0.2;
B=B(B(:,1)==Time,2:3);
B=B(abs(B(:,1))<sqrt(Time/252)*sigma*S0,:)
%}
sigma=0.2;
Exc=0;
Num=2;
B=B(abs(B(:,2))<sqrt(B(:,1)/252)*sigma*S0*Num & B(:,1)>Exc,:);

for i=1:size(B,1)
euro=@(sigma)european_bs(S0,B(i,2)+S0,r,sigma,B(i,1)./252,'put')-B(i,3);
C(i)=fzero(euro,0.5);
end
%scatter3(B(:,1),B(:,2),C);
%vol=fit([B(:,1),B(:,2)],C','poly22');
%vol(100,100)
%plot(vol,[B(:,1),B(:,2)],C')

F = @(var,x)var(1)+var(2)*exp(-var(3)*x(:,1))+var(4)*abs(x(:,2)).*exp(-var(5)*x(:,1))+var(6)*x(:,2).*x(:,2).*exp(-var(7)*x(:,1))+var(8)*abs(x(:,2))+var(9)*x(:,2).*x(:,2)+var(10)*x(:,1)+var(11)*x(:,1).*x(:,1)+var(12)*x(:,1).*exp(-var(13)*x(:,1))+var(14)*x(:,1).*x(:,1).*exp(-var(15)*x(:,1));

data0 = [0.2 0.2 0.004 0.0001 0.004 0.00001 0.004 0.001 0.2 -0.0001 -0.0001 0.01 0.004 0.01 0.004];
opts=optimset('Display','off');
[var] = lsqcurvefit(F,data0,B(:,1:2),C',[],[],opts)
vol=@(x,y)var(1)+var(2)*exp(-var(3)*x)+var(4)*abs(y).*exp(-var(5)*x)+var(6)*y.*y.*exp(-var(7)*x)+var(8)*abs(y)+var(9)*y.*y+var(10)*x+var(11)*x.*x+var(12)*x.*exp(-var(13)*x)+var(14)*x.*x.*exp(-var(15)*x);

scatter3(B(:,1),B(:,2),C,'.');
hold on;
fsurf(vol,[Exc 520 -sqrt(520/252)*sigma*S0*Num sqrt(520/252)*sigma*S0*Num],'ShowContours','on');

%{
hold on;
tri = delaunay(B(:,1),B(:,2)); %plot(B(:,1),B(:,2),'.')
h = trisurf(tri, B(:,1),B(:,2),C'); axis vis3d;
lighting phong;
shading interp;
%}

hold off;
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