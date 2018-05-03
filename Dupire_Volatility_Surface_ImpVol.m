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
SB=B1(B1(:,2)<=MaxK & B1(:,2)>=MinK & B1(:,1)>=MinT & B1(:,1)<=MaxT,:);
dT=10.5/252;
dK=0.05*S01;
Time=MinT:2*dT:MaxT;
Strike=MinK:2*dK:MaxK;
[X,Y] = meshgrid(Time,Strike);


S=scatteredInterpolant(B1(:,1),B1(:,2),B1(:,3),'linear','none');
SgradK=(S(X,Y+dK)-S(X,Y-dK))/(2*dK);
Sgrad2K=(S(X,Y+dK)+S(X,Y-dK)-2*S(X,Y))/(dK^2);
SgradT=(S(X+dT,Y)-S(X,Y))/(dT);
SXY=S(X,Y);

d1=(log(S01./Y)+(r1+0.5*SXY.^2).*X)./(SXY.*sqrt(X));
vol=sqrt((SXY.^2+2*X.*SXY.*SgradT+2*r1*Y.*X.*SXY.*SgradK)./((1+Y.*d1.*sqrt(X).*SgradK).^2+Y.^2.*X.*SXY.*(Sgrad2K-d1.*(SgradK).^2.*sqrt(X))));


V=[];
for i=1:size(vol,1)
   for j=1:size(vol,2)
       if ~isnan(vol(i,j))
           V=[V;[X(i,j),Y(i,j),vol(i,j)]];
       end
   end
end
interp=scatteredInterpolant(V(:,1),V(:,2),V(:,3),'linear','none');
interp2=interp(X,Y);

%%{
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


function euro=european_bs(S0,K,r,sigma0,T,putcall)
d1 = (log(S0/K) + (r + 0.5*sigma0^2)*T)/(sigma0*sqrt(T));
d2 = d1 - sigma0*sqrt(T);
N1 = normcdf(d1);
N2 = normcdf(d2);
if putcall=="call"
    euro = S0*N1 - K*exp(-r*T)*N2;
elseif putcall=="put"
    euro = S0*N1 - K*exp(-r*T)*N2 + K*exp(-r*T) - S0;
end
end