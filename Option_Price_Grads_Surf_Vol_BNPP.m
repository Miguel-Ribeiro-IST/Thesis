clear;

r1 = 0.06;
S01=17099.4;
A = importdata('Data_BNPP.txt','\t',1);
B1=A.data(:,:);

B1(:,2)=B1(:,2)/S01;
S01=1;
B1(:,1)=B1(:,1)/252;
%B1=B1(B1(:,1)<=1 & B1(:,2)<=1.25 & B1(:,2)>=0.74,:);

P1=B1;
for i=1:size(B1,1)
    P1(i,3)=european_bs(S01,B1(i,2),r1,B1(i,3),B1(i,1),"call");
end

%{
for i=1:size(B,1)
euro=@(sigma)european_bs(S0,B(i,2),r,sigma,B(i,1)./252,'put')-B(i,3);
B(i,3)=fzero(euro,0.5);
end
%}


xlabel('time(days)');
ylabel('strike');
zlabel('price');



MinT=min(B1(:,1));
MaxT=max(B1(:,1));
%MinK=min(B1(:,2));
%MaxK=max(B1(:,2));
MinK=0.5;
MaxK=1.70;
SB=P1(P1(:,2)<=MaxK & P1(:,2)>=MinK & P1(:,1)>=MinT & P1(:,1)<=MaxT,:);
dT=42/252;
dK=0.05*S01;
Time=MinT:2*dT:MaxT;
Strike=MinK:2*dK:MaxK;
[X,Y] = meshgrid(Time,Strike);


F=scatteredInterpolant(P1(:,1),P1(:,2),P1(:,3),'linear','none');
FgradK=(F(X,Y+dK)-F(X,Y-dK))/(2*dK);
Fgrad2K=(F(X,Y+dK)+F(X,Y-dK)-2*F(X,Y))/(dK^2);
FgradT=(F(X+dT,Y)-F(X,Y))/(dT);
FS=F(X,Y);

%%{
threshold=0.1;
for i=1:size(Fgrad2K,1)
    for j=1:size(Fgrad2K,2)
        if not(isnan(Fgrad2K(i,j)))
            Fgrad2K(i,j)=max(Fgrad2K(i,j),threshold);
        end
    end
end


%}
%{
for i=1:size(FgradT,1)
    for j=1:size(FgradT,2)
        if not(isnan(FgradT(i,j)))
            FgradT(i,j)=max(FgradT(i,j),-r1*Y(i,j)*FgradK(i,j));
        end
    end
end
%}

vol=2*(FgradT+r1*Y.*FgradK)./(Y.^2.*Fgrad2K);


%%{
ax1 = subplot(2,2,1);
scatter3(ax1,SB(:,1),SB(:,2),SB(:,3),'.');
axis vis3d;
shading interp;
hold on;
surf(ax1,X,Y,FS);
title(ax1,"P");
axis vis3d;
shading interp;
axis([MinT MaxT MinK MaxK])
hold on;
ax2 = subplot(2,2,2);
surf(ax2,X,Y,FgradK);
title(ax2,"dK");
axis vis3d;
shading interp;
axis([MinT MaxT MinK MaxK])
hold on;
ax3 = subplot(2,2,3);
surf(ax3,X,Y,Fgrad2K);
title(ax3,"dK2");
axis vis3d;
shading interp;
axis([MinT MaxT MinK MaxK])
hold on;
ax4 = subplot(2,2,4);
surf(ax4,X,Y,FgradT);
title(ax4,"dT");
axis vis3d;
shading interp;
axis([MinT MaxT MinK MaxK])
%}

%%{
figure
surf(X,Y,vol);
axis vis3d;
shading interp;
%}


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