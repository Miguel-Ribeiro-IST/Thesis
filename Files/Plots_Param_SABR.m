clear;
figure

r=0;
SABRVol=@(K,alpha,rho,nu,beta)sigmaSABR(alpha,rho,nu,beta,K,1,2/12)';
K=(0.4:0.01:1.6);

alpha=0.2;
rho=-0.5;
nu=1.5;
beta=0.5;

%{
alpha1=0.1;
alpha2=0.2;
alpha3=0.3;
for i=1:size(K,2)
    HV1(i)=SABRVol(K(i),alpha1,rho,nu,beta);
    HV2(i)=SABRVol(K(i),alpha2,rho,nu,beta);
    HV3(i)=SABRVol(K(i),alpha3,rho,nu,beta);
end
ttl="\textbf{Dependence on $\alpha$}";
lg={['$\alpha$=',num2str(alpha1),' $\mathrm{yr}^{-1/2}$'],['$\alpha$=',num2str(alpha2),' $\mathrm{yr}^{-1/2}$'],['$\alpha$=',num2str(alpha3),' $\mathrm{yr}^{-1/2}$']};
%}

%{
beta1=0;
beta2=0.5;
beta3=1;
for i=1:size(K,2)
    HV1(i)=SABRVol(K(i),alpha,rho,nu,beta1);
    HV2(i)=SABRVol(K(i),alpha,rho,nu,beta2);
    HV3(i)=SABRVol(K(i),alpha,rho,nu,beta3);
end
ttl="\textbf{Dependence on $\beta$}";
lg={['$\beta$=',num2str(beta1)],['$\beta$=',num2str(beta2)],['$\beta$=',num2str(beta3)]};
%}

%{
rho1=-0.5;
rho2=0;
rho3=0.5;
for i=1:size(K,2)
    HV1(i)=SABRVol(K(i),alpha,rho1,nu,beta);
    HV2(i)=SABRVol(K(i),alpha,rho2,nu,beta);
    HV3(i)=SABRVol(K(i),alpha,rho3,nu,beta);
end
ttl="\textbf{Dependence on $\rho$}";
lg={['$\rho$=',num2str(rho1)],['$\rho$=',num2str(rho2)],['$\rho$=',num2str(rho3)]};
%}


%%{
nu1=0.5;
nu2=1;
nu3=1.5;
for i=1:size(K,2)
    HV1(i)=SABRVol(K(i),alpha,rho,nu1,beta);
    HV2(i)=SABRVol(K(i),alpha,rho,nu2,beta);
    HV3(i)=SABRVol(K(i),alpha,rho,nu3,beta);
end
ttl="\textbf{Dependence on $\nu$}";
lg={['$\nu$=',num2str(nu1),' $\mathrm{yr}^{-1/2}$'],['$\nu$=',num2str(nu2),' $\mathrm{yr}^{-1/2}$'],['$\nu$=',num2str(nu3),' $\mathrm{yr}^{-1/2}$']};
%}

p=plot(K,HV1);
p.LineWidth = 1.5;
hold on;
p=plot(K,HV2,'--');
p.LineWidth = 1.5;
hold on;
p=plot(K,HV3,'-.');
p.LineWidth = 1.5;

xlim([0.4,1.6])
ylim([0,1])
box on;
grid on;
set(gca,'fontsize',12)
xlabel('K/S_0');
ylabel('\sigma_{imp} (yr^{-1/2})')
pbaspect([1.5 1 1])

lgd=legend(lg,'Location','northeast','FontSize',11,'interpreter','latex','FontSize',13);
%title(lgd,ttl,'interpreter','latex','FontSize',13)



function sigma=sigmaSABR(alpha,rho,nu,beta,K,f,T)
z=nu./alpha.*(f.*K).^((1-beta)./2).*log(f./K);
x=log((sqrt(1-2.*rho.*z+z.^2)+z-rho)./(1-rho));
z(z==0,:)=1;
x(x==0,:)=1;

if beta==0
    sigma=alpha.*log(f./K)./(f-K).*(z./x).*(1+T.*(alpha.^2./(24.*f.*K)+(2-3*rho.^2)./24.*nu.^2));
elseif beta==1
    sigma=alpha.*(z./x).*(1+T.*(1/4*rho.*alpha.*nu+1/24*(2-3*rho.^2).*nu.^2));
else
    sigma=alpha./((f.*K).^((1-beta)./2).*(1+(1-beta).^2/24.*(log(f./K)).^2+(1-beta).^4./1920.*(log(f./K)).^4)).*(z./x).*(1+((1-beta).^2./24.*alpha.^2./(f.*K).^(1-beta)+1/4.*rho.*nu.*beta.*alpha./(f.*K).^((1-beta)./2)+(2-3.*rho.^2)./24.*nu.^2).*T);
end
end


function price=european_bs(S0,K,r,sigma,T,putcall)
d1 = (log(S0./K) + (r + 0.5.*sigma.^2).*T)./(sigma.*sqrt(T));
d2 = d1 - sigma.*sqrt(T);
N1 = normcdf(d1);
N2 = normcdf(d2);
if putcall=="call"
    price = S0.*N1 - K.*exp(-r.*T).*N2;
elseif putcall=="put"
    price = S0.*N1 - K.*exp(-r*T).*N2 + K.*exp(-r.*T) - S0;
end
end