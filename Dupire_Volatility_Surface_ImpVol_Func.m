clear;

r1 = 0.06;
S01=17099.4;
A = importdata('Data_BNPP.txt','\t',1);
B1=A.data(:,:);

B1(:,2)=B1(:,2)/S01;
S01=1;
B1(:,1)=B1(:,1)/252;
%B1=B1(B1(:,1)<=0.5,:);

%{
for i=1:size(B,1)
euro=@(sigma)european_bs(S0,B(i,2),r,sigma,B(i,1)./252,'put')-B(i,3);
B(i,3)=fzero(euro,0.5);
end
%}


MinT=min(B1(:,1));
MaxT=max(B1(:,1));
%MinK=min(B1(:,2));
%MaxK=max(B1(:,2));
MinK=0.5;
MaxK=1.70;

dT=10.5/252;
dK=0.05*S01;


func=Dupire(S01,r1,B1,MinT,MaxT,dT,MinK,MaxK,dK);
func(0,[1 1]')
%{
ax1 = subplot(1,2,1);
scatter3(ax1,SB(:,1),SB(:,2),SB(:,3),'.');
axis vis3d;
shading interp;
hold on;
surf(ax1,X,Y,SXY);
title(ax1,"Vol");
axis vis3d;
shading interp;
axis([MinT MaxT MinK MaxK])
hold on;
ax2 = subplot(1,2,2);
surf(ax2,X,Y,interp2);
title(ax2,"vol");
axis vis3d;
shading interp;
axis([MinT MaxT MinK MaxK])
%}


function interp=Dupire(S0,r,B,MinT,MaxT,dT,MinK,MaxK,dK)
B=B(B(:,2)<=MaxK & B(:,2)>=MinK & B(:,1)>=MinT & B(:,1)<=MaxT,:);
Time=MinT:2*dT:MaxT;
Strike=MinK:2*dK:MaxK;
[X,Y] = meshgrid(Time,Strike);
S=scatteredInterpolant(B(:,1),B(:,2),B(:,3),'linear','none');
SgradK=(S(X,Y+dK)-S(X,Y-dK))/(2*dK);
Sgrad2K=(S(X,Y+dK)+S(X,Y-dK)-2*S(X,Y))/(dK^2);
SgradT=(S(X+dT,Y)-S(X,Y))/(dT);
SXY=S(X,Y);

d1=(log(S0./Y)+(r+0.5*SXY.^2).*X)./(SXY.*sqrt(X));
vol=sqrt((SXY.^2+2*X.*SXY.*SgradT+2*r*Y.*X.*SXY.*SgradK)./((1+Y.*d1.*sqrt(X).*SgradK).^2+Y.^2.*X.*SXY.*(Sgrad2K-d1.*(SgradK).^2.*sqrt(X))));


V=[];
for i=1:size(vol,1)
   for j=1:size(vol,2)
       if ~isnan(vol(i,j))
           V=[V;[X(i,j),Y(i,j),vol(i,j)]];
       end
   end
end
interp=scatteredInterpolant(V(:,1),V(:,2),V(:,3),'linear','nearest');
end