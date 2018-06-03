clear;
figure
%{
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
legend({'call','put'},'Location','northeast','FontSize',11)
p.LineWidth=2;
c.LineWidth=2;
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
 yticks(1);
 yticklabels('S_0')
 xticks(1)
  xticklabels('T')
  
xlabel('Time');
ylabel('Stock Price')
%}

%{
M=1000;
S1=ones(1,M);
S2=ones(1,M);
S3=ones(1,M);
r=0.01;
T=1;
dt=T/M;
sigma1=0.1;
sigma2=0.25;
sigma3=0.5;
 for k = 1:M
     Z=randn(1,1);
        S1(:,k+1)=S1(:,k).*(1+r.*dt)+sigma1.*sqrt(dt).*S1(:,k).*Z;
        S2(:,k+1)=S2(:,k).*(1+r.*dt)+sigma2.*sqrt(dt).*S2(:,k).*Z;
        S3(:,k+1)=S3(:,k).*(1+r.*dt)+sigma3.*sqrt(dt).*S3(:,k).*Z;
 end


 hold on;
hold on;
plot(0:dt:1,S2(:,:)','LineWidth',1.25);
hold on;
 plot(0:dt:1,S3(:,:)','LineWidth',1.25);
hold on;

plot(0:dt:1,S1(:,:)','Color',[0. 0.6 0.],'LineWidth',1.25);
hold on;
 yticks(1);
 yticklabels('S_0')
 xticks(1)
  xticklabels('T')
  
xlabel('Time');
ylabel('Stock Price')
%}

%{
smile=@(x)0.02*(x-1).^2+1;
f=fplot(smile,[0,2],'LineWidth',2)
ylim([0.975,1.05])
 xticks(1)
  xticklabels('1.0')
  xlabel('K/S_0');
  yticks([])
ylabel('Implied Volatility')
%}

%{
skew=@(x)0.005*(x-2.3).^2+1;
fplot(skew,[0,2],'LineWidth',2,'Color',[0.4660    0.6740    0.1880])
ylim([0.975,1.05])
 xticks(1)
  xticklabels('1.0')
  xlabel('K/S_0');
  yticks([])
ylabel('Implied Volatility')
%}


%{
M=250;
S1=ones(1,M);
r=0.01;
T=1;
dt=T/M;
sigma1=0.05;
 for k = 1:M
     Z=randn(1,1);
        S1(:,k+1)=S1(:,k).*(1+r.*dt)+sigma1.*sqrt(dt).*S1(:,k).*Z;
 end
 
 delta2=10;
 delta3=50;
 S2=S1(1:delta2:M+1);
 S3=S1(1:delta3:M+1);

plot(0:dt:1,S1(:,:)');
hold on;
plot(0:dt*delta2:1,S2(:,:)');
hold on;
 plot(0:dt*delta3:1,S3(:,:)');
hold on;


 yticks(1);
 yticklabels('S_0')
 xticks(1)
  xticklabels('T')
  
xlabel('Time');
ylabel('Stock Price')
%}