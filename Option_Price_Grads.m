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

B(:,2)=B(:,2)+S0;

for idx=unique(B(:,1))'
    B(B(:,1)==idx,3)=smooth(B(B(:,1)==idx,2),B(B(:,1)==idx,3));
end


%scatter3(B(:,1),B(:,2),C);
%vol=fit([B(:,1),B(:,2)],C','poly22');
%vol(100,100)
%plot(vol,[B(:,1),B(:,2)],C')

scatter3(B(:,1),B(:,2),B(:,3),'.');


hold on;
tri = delaunay(B(:,1),B(:,2)); %plot(B(:,1),B(:,2),'.')
h = trisurf(tri, B(:,1),B(:,2),B(:,3));
axis vis3d;
lighting gouraud;
material dull;
l = light('Position',[0 0 200]);
%l = light('Position',[200 0 0])
shading interp;
%}
xlabel('time(days)');
ylabel('strike');
zlabel('price');

hold on;

T1=9;
K1=190;
dT=5;
dK=1;
F=scatteredInterpolant(B(:,1),B(:,2),B(:,3),'linear','nearest');

price=F(T1,K1);
gradK=(F(T1,K1+dK)-F(T1,K1-dK))/(2*dK);
gradT=(F(T1+dT,K1)-F(T1-dT,K1))/(2*dT);
grad2K=(F(T1,K1+dK)+F(T1,K1-dK)-2*F(T1,K1))/(dK^2);

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

%fsurf(@(x,y)2*vol(x,y),[0 600 -200 200]);
%fsurf(vol,[-200 200 0 600])
%plot(F(var,x,y),[B(:,1),B(:,2)],C')

