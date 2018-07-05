clear;
figure

HestonVol=@(K,T,kappa,nubar,nu0,rho,chi)HestonPrice(K,T,1,0,kappa,nubar,nu0,rho,chi,"vol");
K=(0.4:0.01:1.6);

nubar=0.5^2;
nu0=0.5^2;
rho=0;
chi=1;%eta!!!
kappa=10;

%%{
kappa1=1;
kappa2=25;
kappa3=1000;
for i=1:size(K,2)
    HV1(i)=HestonVol(K(i),2/12,kappa1,nubar,nu0,rho,chi);
    HV2(i)=HestonVol(K(i),2/12,kappa2,nubar,nu0,rho,chi);
    HV3(i)=HestonVol(K(i),2/12,kappa3,nubar,nu0,rho,chi);
end
ttl="\textbf{Dependence on $\kappa$}";
lg={['$\kappa$=',num2str(kappa1)],['$\kappa$=',num2str(kappa2)],['$\kappa$=',num2str(kappa3)]};
%}

%{
nubar1=0.01;
nubar2=0.25;
nubar3=0.5;
for i=1:size(K,2)
    HV1(i)=HestonVol(K(i),2/12,kappa,nubar1,nu0,rho,chi);
    HV2(i)=HestonVol(K(i),2/12,kappa,nubar2,nu0,rho,chi);
    HV3(i)=HestonVol(K(i),2/12,kappa,nubar3,nu0,rho,chi);
end
ttl="\textbf{Dependence on $\overline{\nu}$}";
lg={['$\overline{\nu}$=',num2str(nubar1)],['$\overline{\nu}$=',num2str(nubar2)],['$\overline{\nu}$=',num2str(nubar3)]};
%}

%{
nu01=0.01;
nu02=0.25;
nu03=0.5;
for i=1:size(K,2)
    HV1(i)=HestonVol(K(i),2/12,kappa,nubar,nu01,rho,chi);
    HV2(i)=HestonVol(K(i),2/12,kappa,nubar,nu02,rho,chi);
    HV3(i)=HestonVol(K(i),2/12,kappa,nubar,nu03,rho,chi);
end
ttl="\textbf{Dependence on $\nu_0$}";
lg={['$\nu_0$=',num2str(nu01)],['$\nu_0$=',num2str(nu02)],['$\nu_0$=',num2str(nu03)]};
%}

%{
rho1=-0.5;
rho2=0;
rho3=0.5;
for i=1:size(K,2)
    HV1(i)=HestonVol(K(i),2/12,kappa,nubar,nu0,rho1,chi);
    HV2(i)=HestonVol(K(i),2/12,kappa,nubar,nu0,rho2,chi);
    HV3(i)=HestonVol(K(i),2/12,kappa,nubar,nu0,rho3,chi);
end
ttl="\textbf{Dependence on $\rho$}";
lg={['$\rho$=',num2str(rho1)],['$\rho$=',num2str(rho2)],['$\rho$=',num2str(rho3)]};
%}

%{
chi1=0.5;
chi2=1;
chi3=2;
for i=1:size(K,2)
    HV1(i)=HestonVol(K(i),2/12,kappa,nubar,nu0,rho,chi1);
    HV2(i)=HestonVol(K(i),2/12,kappa,nubar,nu0,rho,chi2);
    HV3(i)=HestonVol(K(i),2/12,kappa,nubar,nu0,rho,chi3);
end
ttl="\textbf{Dependence on $\eta$}";
lg={['$\eta$=',num2str(chi1)],['$\eta$=',num2str(chi2)],['$\eta$=',num2str(chi3)]};
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



function result=HestonPrice(K,T,S0,r,kappa,nubar,nu0,rho,chi,PriceVol)
%Calculate the price of options under Heston using the closed form solution
%We need the characteristic function, depicted in function "CharFuncHeston"
fun1=@(u)real(exp(-1i.*u.*log(K))./(1i.*u.*S0.*exp(r.*T)).*CharFuncHeston(u-1i,T,S0,r,kappa,nubar,nu0,rho,chi));
fun2=@(u)real(exp(-1i.*u.*log(K))./(1i.*u).*CharFuncHeston(u,T,S0,r,kappa,nubar,nu0,rho,chi));
P1=1/2+1/pi.*integral(fun1,0,500);
P2=1/2+1/pi.*integral(fun2,0,500);
call=S0.*P1-exp(-r.*T).*K.*P2; %call option price

if PriceVol=="price"   %if desired output is a price, return this result
    result=call;
else                   %if desired output is a volatility, calculate implied volatility of price
    volatility=@(sigma)european_bs(S0,K,r,sigma,T,"call")-call;
    result=fzero(volatility,1);
end
end

function w=CharFuncHeston(u,T,S0,r,kappa,nubar,nu0,rho,chi)
xi=kappa-chi.*rho.*1i.*u;
d=sqrt(xi.^2+chi.^2.*(u.^2+1i.*u));
A1=(u.^2+1i.*u).*sinh(d.*T/2);
A2=d.*cosh(d.*T./2)+xi.*sinh(d.*T./2);
A=A1./A2;
D=log(d)+(kappa-d).*T./2-log((d+xi)./2+(d-xi)./2.*exp(-d.*T));
w=exp(1i.*u.*(log(S0)+r.*T)-T.*kappa.*nubar.*rho.*1i.*u./chi-nu0.*A+2*kappa.*nubar./chi.^2.*D);
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