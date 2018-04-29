clear;
S0=17099.4;
r=0.06;

A = importdata('Data_BNPP.txt','\t',1);
B=A.data(:,:);


%{
for i=1:size(B,1)
B(i,3)=european_bs(S0,B(i,2),r,B(i,3),B(i,1)./252,'call');
end
%}
B(:,2)=B(:,2)/17099.4;

Month=24;
B1=B(B(:,1)==21*Month,:);




a=0.24;
b=-0.5;
c=0.4;
d=1.1;
e=0.8;
smax=0.9;

%F = @(var,x)min(var(1)+var(2)*(x-var(3)).^2,var(4));
%F =  @(var,x)min(var(1)+(var(2)*(abs(x-var(4))-(x-var(4)))+var(3)*(abs(x-var(4))+(x-var(4)))).*(x-var(4)).^2,e);
%F =  @(var,x)min(var(1)+var(2)*exp(-var(3)*(x-var(4)))+var(5)*(x-var(4)),smax);
F =  @(var,x)min(var(1)+(-var(2)*(abs(x-var(4))-(x-var(4)))+var(3)*(abs(x-var(4))+(x-var(4)))).*(x-var(4)),smax);
%F =  @(var,x)min(var(1)+(var(2)*(abs(x-var(4))-(x-var(4)))+var(3)*(abs(x-var(4))+(x-var(4)))),smax);

data0 = [a b c d e];
data0 = [a b d e];
opts=optimset('Display','off');
[var] = lsqcurvefit(F,data0,B1(:,2),B1(:,3),[],[],opts)
%vol=@(x)min(var(1)+(var(2)*(abs(x-var(4))-(x-var(4)))+var(3)*(abs(x-var(4))+(x-var(4)))).*(x-var(4)).^2,e);
vol=@(x)min(var(1)+(-var(2)*(abs(x-var(4))-(x-var(4)))+var(3)*(abs(x-var(4))+(x-var(4)))).*(x-var(4)),smax);
%vol=@(x)min(var(1)+var(2)*(x-var(3)).^2,var(4));

var2= [a b c d e];
vol2=@(x)(0.2*(abs(x-var2(4))-(x-var2(4)))+0.3*(abs(x-var2(4))+(x-var2(4))));


scatter(B1(:,2),B1(:,3));
hold on;
%{
scatter(B2(:,2),B2(:,3));
hold on;

scatter(B3(:,2),B3(:,3));
hold on;
scatter(B6(:,2),B6(:,3));
hold on;
%}
fplot(vol,[0.4,1.7]);
hold on;
fplot(vol2, [0.4,1.7]);
hold on;

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