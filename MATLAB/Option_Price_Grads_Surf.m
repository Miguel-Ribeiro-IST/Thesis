clear;
S0=177.12;
r=0.06;

A = importdata('Options_AAPL.txt',';',1);
B=A.data(:,:);

sigma=0.2;
Exc=10;
Excm=200;
Num=1;
%B=B(abs(B(:,2))<sqrt(B(:,1)/252)*sigma*S0*Num & B(:,1)>Exc & B(:,1)<Excm,:);


B(:,2)=B(:,2)+S0;


for idx=unique(B(:,1))'
    B(B(:,1)==idx,3)=smooth(B(B(:,1)==idx,2),B(B(:,1)==idx,3));
end


xlabel('time(days)');
ylabel('strike');
zlabel('price');



dT=1;
dK=1;
Time=0:2*dT:(max(B(:,1)));
TimeP=Time+dT;
TimeM=Time-dT;
Strike=0:2*dK:(max(B(:,2)));
StrikeP=Strike+dK;
StrikeM=Strike-dK;
[X,Y] = meshgrid(Time,Strike);
[XKP,YKP] = meshgrid(Time,StrikeP);
[XKM,YKM] = meshgrid(Time,StrikeM);
[XTP,YTP] = meshgrid(TimeP,Strike);
[XTM,YTM] = meshgrid(TimeM,Strike);


F=scatteredInterpolant(B(:,1),B(:,2),B(:,3),'linear','none');
FgradK=(F(XKP,YKP)-F(XKM,YKM))/(2*dK);
Fgrad2K=(F(XKP,YKP)+F(XKM,YKM)-2*F(X,Y))/(dK^2);
FgradT=(F(XTP,YTP)-F(XTM,YTM))/(2*dT);
FS=F(X,Y);

vol=(2*FgradT+r*Y.*FgradK)./(Y.^2.*Fgrad2K);

threshold=0.3;

for i=1:size(vol,1)
    for j=1:size(vol,2)
        if not(isnan(vol(i,j)))
        vol(i,j)=min(vol(i,j),threshold);
        end
    end
end

%scatter3(B(:,1),B(:,2),B(:,3),'.');
%hold on;
s=surf(X,Y,FS);
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

