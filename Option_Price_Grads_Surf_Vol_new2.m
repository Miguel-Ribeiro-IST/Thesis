clear;
S0=167.65;
r=0.06;

A = importdata('AAPL.txt','\t',1);
B=A.data(:,:);


sigma=0.2;

%%{
for idx=unique(B(:,1))' 
C=zeros(size(B(B(:,1)==idx,3),1),1);
C(:,1)=B(B(:,1)==idx,3);
D=zeros(size(B(B(:,1)==idx,3),1),1);
D(:,1)=B(B(:,1)==idx,3);
U=zeros(size(B(B(:,1)==idx,3),1),1);
U(:,1)=B(B(:,1)==idx,3);
while 1
    bk=0;
for i=3:(size(B(B(:,1)==idx,3),1)-1)
    if C(i,1)>C(i-1,1)
        U(i-1,1)=(U(i,1)+U(i-2,1))/2;
        D(i,1)=(D(i-1,1)+D(i+1,1))/2;
        C([i-1 i])=C([i i-1]);
        bk=1;
    end
end
if bk==0
    break
end
end
B(B(:,1)==idx,3)=(U(:,1)+D(:,1))/2;
end
%}

%%{
for i=1:size(B,1)
euro=@(sigma)european_bs(S0,B(i,2),r,sigma,B(i,1)./365,'call')-B(i,3);
B(i,3)=fzero(euro,0.5);
end
%}


%{
for idx=unique(B(:,1))'
    B(B(:,1)==idx,3)=smooth(B(B(:,1)==idx,2),B(B(:,1)==idx,3));
end
%}

xlabel('time(days)');
ylabel('strike');
zlabel('price');


MinT=20;
MaxT=100;
MinK=157.4;
MaxK=177.6;
SB=B(B(:,2)<MaxK & B(:,2)>MinK & B(:,1)>MinT & B(:,1)<MaxT,:);
dT=5;
dK=0.5;
Time=MinT:2*dT:MaxT;
Strike=MinK:2*dK:MaxK;
[X,Y] = meshgrid(Time,Strike);


F=scatteredInterpolant(B(:,1),B(:,2),B(:,3),'natural','none');
FgradK=(F(X,Y+dK)-F(X,Y-dK))/(2*dK);
Fgrad2K=(F(X,Y+dK)+F(X,Y-dK)-2*F(X,Y))/(dK^2);
FgradT=(F(X+dT,Y)-F(X-dT,Y))/(2*dT);
FS=F(X,Y);


vol=(2*FgradT+r*Y.*FgradK)./(Y.^2.*Fgrad2K);

%{
for i=1:size(Fgrad2K,1)
    for j=1:size(Fgrad2K,2)
        if not(isnan(Fgrad2K(i,j)))
        Fgrad2K(i,j)=min(Fgrad2K(i,j),0);
        end
    end
end
%}

%%{
scatter3(SB(:,1),SB(:,2),SB(:,3),'.');
hold on;
s=surf(X,Y,FS);
%}
%s=surf(X,Y,vol);
%s=surf(X,Y,FgradK);
%s=surf(X,Y,Fgrad2K);
%s=surf(X,Y,FgradT);
%axis([0 max(B(:,1)) 0 max(B(:,2)) 0 0.5]);
axis vis3d;
shading interp;
%{
price=F(T1,K1);
gradK=(F(T1,K1+dK)-F(T1,K1-dK))/(2*dK);
gradT=(F(T1+dT,K1)-F(T1-dT,K1))/(2*dT);
grad2K=(F(T1,K1+dK)+F(T1,K1-dK)-2*F(T1,K1))/(dK^2)
%}

%{
scatter3(T1,K1,price,'.','red');
hold on;
syms t
xT = t/gradT+T1;
yT = t*0+K1;
zT = t+price;
fplot3(xT,yT,zT,'red',[-price max(B(:,3))-price])

xK = t*0+T1;
yK = t/gradK+K1;
zK = t+price;
fplot3(xK,yK,zK,'red',[-price max(B(:,3))-price])
axis([0 max(B(:,1)) 0 max(B(:,2)) 0 max(B(:,3))])

hold off;

loc_vol=sqrt((2*gradT+K1*r*gradK)/(K1^2*grad2K))
%}

%fsurf(@(x,y)2*vol(x,y),[0 600 -200 200]);
%fsurf(vol,[-200 200 0 600])
%plot(F(var,x,y),[B(:,1),B(:,2)],C')

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