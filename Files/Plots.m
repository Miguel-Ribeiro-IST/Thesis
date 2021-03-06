clear;
%figure
%%{
figure
K=1;
syms S;
call=piecewise(S<K,0,S>K,S-K);
put=piecewise(S>K,0,S<K,K-S);
c=fplot(call,[0,2]);
hold on;
p=fplot(put,[0,2],'--');
xlabel('S(T)');
ylabel('Payoff')
yticks([]);
xticks([1]);
xticklabels({'K'})
legend({'Call','Put'},'Location','northeast','FontSize',11)
p.LineWidth=2;
c.LineWidth=2;
pbaspect([1.5 1 1])
set(gca,'fontsize',12)
%}

%{
M=1000;
N=3;
S=ones(N,M);
r=0.01;
T=1;
dt=T/M;
sigma=0.1;
 for k = 1:M
     
        S(:,k+1)=S(:,k).*(1+r.*dt)+sigma.*sqrt(dt).*S(:,k).*randn(N,1);
 end

p=plot(0:dt:1,S(:,:)','LineWidth',1.25);
% yticks(1);
% yticklabels('S_0')
% xticks(1)
%  xticklabels('T')
ylim([0.85,1.15])
xlabel('Time(yr)');
ylabel('Stock Price(�)')
pbaspect([1.5 1 1])
set(gca,'fontsize',12)

%}

%{
M=1000;
S1=ones(1,M);
S2=ones(1,M);
S3=ones(1,M);
r=0.01;
T=1;
dt=T/M;
sigma1=0.05;
sigma2=0.1;
sigma3=0.2;
 for k = 1:M
     Z=randn(1,1);
        S1(:,k+1)=S1(:,k).*(1+r.*dt)+sigma1.*sqrt(dt).*S1(:,k).*Z;
        S2(:,k+1)=S2(:,k).*(1+r.*dt)+sigma2.*sqrt(dt).*S2(:,k).*Z;
        S3(:,k+1)=S3(:,k).*(1+r.*dt)+sigma3.*sqrt(dt).*S3(:,k).*Z;
 end



plot(0:dt:1,S2(:,:)','-.','LineWidth',1.25);
hold on;
 plot(0:dt:1,S3(:,:)','LineWidth',1.25);
hold on;

plot(0:dt:1,S1(:,:)',':','LineWidth',1.25);

h = get(gca,'Children');
legend([h(1) h(3) h(2)],{'\sigma=0.05 yr^{-0.5}','\sigma=0.1 yr^{-0.5}','\sigma=0.2 yr^{-0.5}'},'Location','northeast','FontSize',11)
  

set(gca,'fontsize',12)

ylim([0.75,1.25])
xlabel('Time(yr)');
ylabel('Stock Price(�)')
pbaspect([1.5 1 1])
%}

%{
smile=@(x)0.02*(x-1).^2+1;
fplot(smile,[0,2],'LineWidth',2)
ylim([0.975,1.05])
 xticks(1)
  xticklabels('1.0')
  xlabel('K/S_0');
  yticks([])
ylabel('Implied Volatility')

set(gca,'fontsize',12)
pbaspect([1.5 1 1])
%}

%{
skew=@(x)0.005*(x-2.3).^2+1;
fplot(skew,[0,2],'LineWidth',2,'Color',[0.8500    0.3250    0.0980])
ylim([0.975,1.05])
 xticks(1)
  xticklabels('1.0')
  xlabel('K/S_0');
  yticks([])
ylabel('Implied Volatility')

set(gca,'fontsize',12)

pbaspect([1.5 1 1])
%}


%{

M=500;
S1=ones(1,M);
r=0.01;
T=1;
dt=T/M;
sigma1=0.2;
 for k = 1:M
     Z=randn(1,1);
        S1(:,k+1)=S1(:,k).*(1+r.*dt)+sigma1.*sqrt(dt).*S1(:,k).*Z;
 end
 
 delta2=10;
 delta3=50;
 S2=S1(1:delta2:M+1);
 S3=S1(1:delta3:M+1);

plot(0:dt:1,S1(:,:)','LineWidth',1.25);
hold on;
plot(0:dt*delta2:1,S2(:,:)','LineWidth',1.5);
hold on;
 plot(0:dt*delta3:1,S3(:,:)','LineWidth',2);
hold on;
legend({'\Delta t=0.002 yr','\Delta t=0.02 yr','\Delta t=0.1 yr'},'FontSize',11)


set(gca,'fontsize',12)

ylim([0.75,1.25])
xlabel('Time(yr)');
ylabel('Stock Price(�)')
pbaspect([1.5 1 1])

%}

%{
Val=[72, 80, 81, 88, 94, 96, 100, 111, 128, 142, 169, 197, 220, 259, 283, ...
299, 373, 418, 507, 586, 673, 598, 595, 603, 582, 601, 707, 648, 641, ...
635, 696, 710, 691, 628, 551, 493, 553, 482, 542];
startDate = datenum('06-01-1998');
endDate = datenum('06-01-2017');
xData = linspace(startDate,endDate,39);
plot(xData,Val,'LineWidth',2)
datetick('x','mmm-yy')
xlim(1.0e+05*[7.29907  7.36847])
xticks(xData(1:5:end))
grid on;
ylabel('OTC Market Size (Trillion $)')
xlabel('Date')
set(gca,'fontsize',12)

pbaspect([1.5 1 1])
%}

%{
A = importdata('Data_BNPP.txt','\t',1);
B=A.data(:,:);

%%%%%%%%%%%%%%%%%%%%  INPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%
S0=17099.4;        %initial stock price
r = 0;
matur=4;
%%%%%%%%%%%%%      ORIGINAL DATA MODIFICATIONS     %%%%%%%%%%%%
B(:,2)=B(:,2)/S0;     %normalize strike prices
S0=1;
B(:,1)=B(:,1)/252;    %convert maturities from days to years

times=unique(B(:,1));
T = times(matur);

C=B(B(:,1)==T,2:3);


scatter(C(:,1),C(:,2),100,'x','LineWidth',1.5);
    xlim([0.4,1.6])
    ylim([0,1])
    box on;
set(gca,'fontsize',12)

xlabel('K/S_0');
ylabel('\sigma_{imp}(yr^{-1/2})')
grid on;
pbaspect([1.5 1 1])
%}


%{
A = importdata('Data_BNPP.txt','\t',1);
B=A.data(:,:);

%%%%%%%%%%%%%%%%%%%%  INPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%
S0=17099.4;        %initial stock price
r = 0;
matur=6;
%%%%%%%%%%%%%      ORIGINAL DATA MODIFICATIONS     %%%%%%%%%%%%
B(:,2)=B(:,2)/S0;     %normalize strike prices
S0=1;
B(:,1)=B(:,1)/252;    %convert maturities from days to years

times=unique(B(:,1));
T = times(matur);
format long
P=B;
for i=1:size(B,1)
    P(i,3)=european_bs(S0,B(i,2),r,B(i,3),B(i,1),"call");
end

C=P(P(:,1)==T,2:3);


scatter(C(:,1),C(:,2),100,'x','LineWidth',1.5);
    xlim([0.4,1.6])
    ylim([0,0.6])
    box on;
set(gca,'fontsize',12)

xlabel('K/S_0');
ylabel('Option Price (�)')
grid on;
pbaspect([1.5 1 1])

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
%}

%{
weight=@(x)(1-abs(x-1)).^2;
fplot(weight,[0,2],'LineWidth',2)
x=[0.5,0.75,0.9,1,1.1,1.25,1.5];
y=weight(x)'
ylim([0,1.2])
  xlabel('K/S_0');
ylabel('Weight')

set(gca,'fontsize',12)
grid on;
pbaspect([1.5 1 1])
%}

%{
figure
format long
T=1/12;
r=0.01;
S0=1;
m=0.9;
d=1000;

K1=0.9;
P=linspace(S0-K1*exp(-r*T),0.3,d);
vol=zeros(1,size(P,2));
for i=1:size(P,2)
volatility=@(sigma)european_bs(S0,K1,r,sigma,T,"call")-P(i);
vol(i)=fzero(volatility,0.000001);
end
plot(P,vol,'LineWidth',2,'Color',[0.9500    0.2250    0.0580]);
hold on;

K2=0.75;
P=linspace(S0-K2*exp(-r*T),0.4,d);
vol=zeros(1,size(P,2));
for i=1:size(P,2)
volatility=@(sigma)european_bs(S0,K2,r,sigma,T,"call")-P(i);
vol(i)=fzero(volatility,0.000001);
end
plot(P,vol,'-.','LineWidth',2,'Color',[0.0010    0.60    0.8330]);
hold on;

K3=0.5;
P=linspace(S0-K3*exp(-r*T),0.55,d);
vol=zeros(1,size(P,2));
for i=1:size(P,2)
volatility=@(sigma)european_bs(S0,K3,r,sigma,T,"call")-P(i);
vol(i)=fzero(volatility,0.000001);
end
plot(P,vol,':','LineWidth',2.25,'Color',[0.4660    0.6740    0.1880]);

xlim([0,0.6])
ylim([0,2])
%set(gca, 'YScale', 'log')

xlabel('C_{mkt}(�)');
ylabel('\sigma_{imp}(yr^{-1/2})')
%xticks(S0-K*exp(-r*T));
%xticklabels({'S_0-Ke^{(-rT)}'})
pbaspect([1.5 1 1])
set(gca,'fontsize',12)
grid on;
legend({['K=',num2str(K1),'0 �'],['K=',num2str(K2),' �'],['K=',num2str(K3),'0 �']},'Location','northwest','FontSize',11)


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
%}

%{
figure
format long
T=1/12;
r=0.01;
S0=1;
m=0.1;
d=1000;

K1=1.1;
P=linspace(0,m,d);
vol=zeros(1,size(P,2));
for i=1:size(P,2)
volatility=@(sigma)european_bs(S0,K1,r,sigma,T,"call")-P(i);
vol(i)=fzero(volatility,0.000001);
end
plot(vol,P,'LineWidth',2,'Color',[0.9500    0.2250    0.0580]);
hold on;

K2=1.25;
vol=zeros(1,size(P,2));
for i=1:size(P,2)
volatility=@(sigma)european_bs(S0,K2,r,sigma,T,"call")-P(i);
vol(i)=fzero(volatility,0.000001);
end
plot(vol,P,'-.','LineWidth',2,'Color',[0.0010    0.60    0.8330]);
hold on;

K3=1.5;
vol=zeros(1,size(P,2));
for i=1:size(P,2)
volatility=@(sigma)european_bs(S0,K3,r,sigma,T,"call")-P(i);
vol(i)=fzero(volatility,0.000001);
end
plot(vol,P,':','LineWidth',2.25,'Color',[0.4660    0.6740    0.1880]);

xlim([0,1])
ylim([0,m])
%set(gca, 'YScale', 'log')

ylabel('C_{mkt}(�)');
xlabel('\sigma_{imp}(yr^{-1/2})')
%xticks(S0-K*exp(-r*T));
%xticklabels({'S_0-Ke^{(-rT)}'})
pbaspect([1.5 1 1])
set(gca,'fontsize',12)
grid on;
legend({['K=',num2str(K1),'0 �'],['K=',num2str(K2),' �'],['K=',num2str(K3),'0 �']},'Location','northwest','FontSize',11)


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
%}

%{
figure
format long
r=0.01;
S0=1;
d=1000;
sigma=0.3;

K=linspace(0.4,1.6,10000);
veg=zeros(1,size(K,2));
pric=zeros(1,size(K,2));

T=3/12;
for i=1:size(K,2)
pric(i)=european_bs(S0,K(i),r,sigma,T,"call");
veg(i)=vega(S0,K(i),r,sigma,T)/pric(i)*sigma;
end
plot(K,veg,'LineWidth',2);
hold on;

T=2/12;
for i=1:size(K,2)
pric(i)=european_bs(S0,K(i),r,sigma,T,"call");
veg(i)=vega(S0,K(i),r,sigma,T)/pric(i)*sigma;
end
plot(K,veg,'-.','LineWidth',2);
hold on;

T=1/12;
for i=1:size(K,2)
pric(i)=european_bs(S0,K(i),r,sigma,T,"call");
veg(i)=vega(S0,K(i),r,sigma,T)/pric(i)*sigma;
end
plot(K,veg,':','LineWidth',2);
hold on;

xlim([0.4,1.6])
%ylim([0,0.15])
%set(gca, 'YScale', 'log')

xlabel('K/S_0');
ylabel('Relative Change')
%xticks(S0-K*exp(-r*T));
%xticklabels({'S_0-Ke^{(-rT)}'})
pbaspect([1.5 1 1])
set(gca,'fontsize',12)
grid on;
legend({'T=63 days','T=42 days','T=21 days'},'Location','northwest','FontSize',11)


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

function res=vega(S0,K,r,sigma,T)
d1 = (log(S0./K) + (r + 0.5.*sigma.^2).*T)./(sigma.*sqrt(T));
N=@(x)1/(sqrt(2*pi))*exp(-x.^2/2);
res=S0*sqrt(T)*N(d1);
end
%}

%{
figure
format long
r=0.01;
S0=1;
d=1000;
sigma=0.3;

K=linspace(0.4,1.6,100000);
veg=zeros(1,size(K,2));
pric=zeros(1,size(K,2));

T=3/12;
for i=1:size(K,2)
pric(i)=european_bs(S0,K(i),r,sigma,T,"call");
veg(i)=vega(S0,K(i),r,sigma,T);
end
plot(K,veg,'LineWidth',2);
hold on;

T=2/12;
for i=1:size(K,2)
pric(i)=european_bs(S0,K(i),r,sigma,T,"call");
veg(i)=vega(S0,K(i),r,sigma,T);
end
plot(K,veg,'-.','LineWidth',2);
hold on;

T=1/12;
for i=1:size(K,2)
pric(i)=european_bs(S0,K(i),r,sigma,T,"call");
veg(i)=vega(S0,K(i),r,sigma,T);
end
plot(K,veg,':','LineWidth',2);
hold on;

xlim([0.4,1.6])
%ylim([0,0.15])
%set(gca, 'YScale', 'log')

xlabel('K/S_0');
ylabel('Vega (�.yr^{1/2})')
%xticks(S0-K*exp(-r*T));
%xticklabels({'S_0-Ke^{(-rT)}'})
pbaspect([1.5 1 1])
set(gca,'fontsize',12)
grid on;
legend({'T=63 days','T=42 days','T=21 days'},'Location','northeast','FontSize',11)


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

function res=vega(S0,K,r,sigma,T)
d1 = (log(S0./K) + (r + 0.5.*sigma.^2).*T)./(sigma.*sqrt(T));
N=@(x)1/(sqrt(2*pi))*exp(-x.^2/2);
res=S0*sqrt(T)*N(d1);
end
%}



%%{
figure
format long
r=0.01;
S0=1;
d=1000;
sigma=0.3;

K=linspace(0.4,1.6,1000);
pric=zeros(1,size(K,2));
pric2=zeros(1,size(K,2));

T=6/12;
for i=1:size(K,2)
pric=european_bs(S0,K,r,sigma,T,"call");
end
plot(K,pric,'LineWidth',2);
hold on;


for i=1:size(K,2)
pric2=european_bs(S0,K,r,sigma,T,"put");
end
plot(K,pric2,'LineWidth',2);
hold on;

T=0;
for i=1:size(K,2)
pric=european_bs(S0,K,r,sigma,T,"call");
end
plot(K,pric,'--','LineWidth',2);
hold on;


for i=1:size(K,2)
pric2=european_bs(S0,K,r,sigma,T,"put");
end
plot(K,pric2,'--','LineWidth',2);
hold on;



xlim([0.4,1.6])
ylim([0,0.6])

xlabel('K/S_0');
ylabel('Value (�)')
pbaspect([1.5 1 1])
set(gca,'fontsize',12)
grid on;
legend({'Call (at Inception)','Put (at Inception)','Call (at Maturity)','Put (at Maturity)'},'Location','northeast','FontSize',11)


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

%}
