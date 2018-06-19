%clear;

r = 0.00;
S0=17099.4;
A = importdata('Data_BNPP.txt','\t',1);
B=A.data(:,:);
matur=4;

B(:,2)=B(:,2)/S0;
S0=1;
B(:,1)=B(:,1)/252;
times=unique(B(:,1));
B=B(B(:,1)<=times(matur),:);


MinT=times(1);
MaxT=times(matur);
MinK=0.4;
MaxK=1.6;
SB=B(B(:,2)<=MaxK & B(:,2)>=MinK & B(:,1)>=MinT & B(:,1)<=MaxT,:);
dT=10.5/252;
dK=0.05*S0;
Time=MinT:2*dT:MaxT;
Strike=MinK:2*dK:MaxK;
[X,Y] = meshgrid(Strike,Time);


S=scatteredInterpolant(B(:,2),B(:,1),B(:,3),'linear','linear');
SgradK=(S(X+dK,Y)-S(X-dK,Y))/(2*dK);
Sgrad2K=(S(X+dK,Y)+S(X-dK,Y)-2*S(X,Y))/(dK^2);
SgradT=(S(X,Y+dT)-S(X,Y))/(dT);
SXY=S(X,Y);

d1=(log(S0./X)+(r+0.5*SXY.^2).*Y)./(SXY.*sqrt(Y));
vol=sqrt((SXY.^2+2*Y.*SXY.*SgradT+2*r*X.*Y.*SXY.*SgradK)./((1+X.*d1.*sqrt(Y).*SgradK).^2+X.^2.*Y.*SXY.*(Sgrad2K-d1.*(SgradK).^2.*sqrt(Y))));


V=[];
for i=1:size(vol,1)
    for j=1:size(vol,2)
        if ~isnan(vol(i,j))
            V=[V;[X(i,j),Y(i,j),vol(i,j)]];
        end
    end
end
interp=scatteredInterpolant(V(:,1),V(:,2),V(:,3),'linear','linear');
interp2=interp(X,Y);

%%{
figure
scatter3(SB(:,2),SB(:,1),SB(:,3),30,'LineWidth',0.5,'MarkerEdgeColor','k','MarkerFaceColor',[0.3010    0.7450    0.9330]);
hold on;
surf(X,Y,SXY);
s.EdgeAlpha=0.6;
s.FaceAlpha=0.85;
shading interp
axis([MinK MaxK MinT MaxT 0 1])
box on;
grid on;
xlabel('K/S_0');
ylabel('T (days)');
zlabel('\sigma_{imp} (yr^{-1/2})')
yticks([1/12,2/12,3/12,4/12,5/12,6/12])
yticklabels({'21','42','63','84','105','126'})
pbaspect([1 1.5 1])
set(gca,'fontsize',11)
view(40,35)
M = view(gca);
R = M(1:3,1:3);
x = R*[1;0;0];
y = R*[0;1;0];
z = R*[0;0;1];
set(get(gca,'XLabel'),'rotation',360/(2*pi)*atan(x(2)/x(1)))
set(get(gca,'YLabel'),'rotation',360/(2*pi)*atan(y(2)/y(1)))

figure
contourf(X,Y,SXY,25)
pbaspect([1.5 1 1])
xlim([0.4,1.6])
ylim([1/12,0.5])
xlabel('K/S_0');
ylabel('T (days)');
yticks([1/12,2/12,3/12,4/12,5/12,6/12])
yticklabels({'21','42','63','84','105','126'})
box on;
colorbar;
set(gca,'fontsize',12)




figure
surf(X,Y,interp2);
axis([MinK MaxK MinT MaxT 0 1])
s.EdgeAlpha=0.6;
s.FaceAlpha=0.85;
shading interp
axis([MinK MaxK MinT MaxT])
box on;
grid on;
xlabel('K/S_0');
ylabel('T (days)');
zlabel('\sigma_{loc} (yr^{-1/2})')
yticks([1/12,2/12,3/12,4/12,5/12,6/12])
yticklabels({'21','42','63','84','105','126'})
pbaspect([1 1.5 1])
set(gca,'fontsize',11)
view(40,35)
M = view(gca);
R = M(1:3,1:3);
x = R*[1;0;0];
y = R*[0;1;0];
z = R*[0;0;1];
set(get(gca,'XLabel'),'rotation',360/(2*pi)*atan(x(2)/x(1)))
set(get(gca,'YLabel'),'rotation',360/(2*pi)*atan(y(2)/y(1)))

figure
contourf(X,Y,interp2,25)
pbaspect([1.5 1 1])
xlabel('K/S_0');
ylabel('T (days)');
xlim([0.4,1.6])
ylim([1/12,0.5])
yticks([1/12,2/12,3/12,4/12,5/12,6/12])
yticklabels({'21','42','63','84','105','126'})
box on;
colorbar;
set(gca,'fontsize',12)


function euro=european_bs(S0,K,r,sigma0,T,putcall)
d1 = (log(S0/K) + (r + 0.5*sigma0^2)*T)./(sigma0*sqrt(T));
d2 = d1 - sigma0*sqrt(T);
N1 = normcdf(d1);
N2 = normcdf(d2);
if putcall=="call"
    euro = S0*N1 - K*exp(-r*T)*N2;
elseif putcall=="put"
    euro = S0*N1 - K*exp(-r*T)*N2 + K*exp(-r*T) - S0;
end
end