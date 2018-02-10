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
Excm=2000;
Num=20;
B=B(abs(B(:,2))<sqrt(B(:,1)/252)*sigma*S0*Num & B(:,1)>Exc & B(:,1)<Excm,:);

B(:,2)=B(:,2)+S0;

%{
for idx=unique(B(:,1))'
    B(B(:,1)==idx,3)=smooth(B(B(:,1)==idx,2),B(B(:,1)==idx,3));
end
%}

%scatter3(B(:,1),B(:,2),C);
%vol=fit([B(:,1),B(:,2)],C','poly22');
%vol(100,100)
%plot(vol,[B(:,1),B(:,2)],C')
gradK=zeros(1,3);
for time=unique(B(:,1))'
    dK1=diff(B(B(:,1)==time,2));
    dK2=dK1;
    dK1=cat(1,0,dK1);
    dK2=cat(1,dK2,0);
    dK=dK1+dK2;
    dK=dK(2:end-1);
    dP1=diff(B(B(:,1)==time,3));
    dP2=dP1;
    dP1=cat(1,0,dP1);
    dP2=cat(1,dP2,0);
    dP=dP1+dP2;
    dP=dP(2:end-1);
    tmp=B(B(:,1)==time,1:2);
    tmp=tmp(2:end-1,:);
    gradK=cat(1,gradK,[tmp,dP./dK]);
end
gradK=gradK(2:end,:);


grad2K=zeros(1,3);
for time=unique(B(:,1))'
    dK1=diff(B(B(:,1)==time,2));
    dK2=dK1;
    dK1=cat(1,0,dK1);
    dK2=cat(1,dK2,0);
    dK=dK1+dK2;
    dK=dK(2:end-1);
    dP1=diff(B(B(:,1)==time,3));
    dP2=dP1;
    dP1=cat(1,0,dP1);
    dP2=cat(1,dP2,0);
    dP=dP1-dP2;
    dP=dP(2:end-1);
    tmp=B(B(:,1)==time,1:2);
    tmp=tmp(2:end-1,:);
    grad2K=cat(1,grad2K,[tmp,dP./(dK.^2)]);
end
grad2K=grad2K(2:end,:);


for idx=unique(gradK(:,1))'
    gradK(gradK(:,1)==idx,3)=smooth(gradK(gradK(:,1)==idx,2),gradK(gradK(:,1)==idx,3));
end
%}

for idx=unique(grad2K(:,1))'
    grad2K(grad2K(:,1)==idx,3)=smooth(grad2K(grad2K(:,1)==idx,2),grad2K(gradK(:,1)==idx,3));
end
%}

%{
hold on;
tri = delaunay(gradK(:,1),gradK(:,2)); %plot(B(:,1),B(:,2),'.')
h = trisurf(tri, gradK(:,1),gradK(:,2),gradK(:,3)); axis vis3d;
lighting gouraud;
material dull;
l = light('Position',[200 0 1]);
%l = light('Position',[200 0 0])
shading interp;
%}

Time=0:5:(max(B(:,1)));
Strike=0:1:(max(B(:,2)));
[X,Y] = meshgrid(Time,Strike);

Fk=scatteredInterpolant(gradK(:,1),gradK(:,2),gradK(:,3),'natural','none');
F2k=scatteredInterpolant(grad2K(:,1),grad2K(:,2),grad2K(:,3),'natural','none');

gradKm=Fk(X,Y);
grad2Km=F2k(X,Y);

%scatter3(gradK(:,1),gradK(:,2),gradK(:,3),'.');
scatter3(grad2K(:,1),grad2K(:,2),grad2K(:,3),'.');
s=meshc(X,Y,grad2Km);
axis([0 max(B(:,1)) 0 max(B(:,2)) -0.005 0.05])
%set(s,'LineStyle','none')




%{
hold on;
tri = delaunay(grad2K(:,1),grad2K(:,2)); %plot(B(:,1),B(:,2),'.')
h = trisurf(tri, grad2K(:,1),grad2K(:,2),grad2K(:,3)); axis vis3d;
lighting gouraud;
material dull;
l = light('Position',[200 0 1]);
%l = light('Position',[200 0 0])
shading interp;
%}


%{
dT=5;
dK=0.2;
marginT=10;
marginK=10;
Time=0:5:(max(B(:,1)));
Strike=0:1:(max(B(:,2)));
[X,Y] = meshgrid(Time,Strike);



%F=scatteredInterpolant(B(:,1),B(:,2),B(:,3),'linear','nearest');
F=scatteredInterpolant(B(:,1),B(:,2),B(:,3),'natural','none');

gradK=(F(X,Y+dK)-F(X,Y-dK)+(F(X+marginT,Y)+F(X-marginT,Y)+F(X,Y+marginK)+F(X,Y-marginK))*0)/(2*dK);
gradT=(F(X+dT,Y)-F(X-dT,Y))/(2*dT);
gradK2=(F(X,Y+dK)+F(X,Y-dK)-2*F(X,Y))/(dK^2);
vols=(2*gradT+r*Y.*gradK)./(Y.^2.*gradK2);
    %gradK=(F(T1,K1+dK)-F(T1,K1-dK))/(2*dK)
%gradT=(F(T1+dT,K1)-F(T1-dT,K1))/(2*dT)


%{
tri = delaunay(B(:,1),B(:,2)); %plot(B(:,1),B(:,2),'.')
h = trisurf(tri, B(:,1),B(:,2),B(:,3)); axis vis3d;
lighting gouraud;
material dull;
l = light('Position',[200 0 1]);
%l = light('Position',[200 0 0])
shading interp;
%}

s=meshc(X,Y,vols);
axis vis3d;
%set(s,'LineStyle','none')
%axis([0 max(B(:,1)) 0 max(B(:,2)) 0 1])
%{
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
%}

%fsurf(@(x,y)2*vol(x,y),[0 600 -200 200]);
%fsurf(vol,[-200 200 0 600])
%plot(F(var,x,y),[B(:,1),B(:,2)],C')

%}