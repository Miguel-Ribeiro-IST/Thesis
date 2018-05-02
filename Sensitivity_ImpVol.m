clear;

r1 = 0.06;
S01=17099.4;
A = importdata('Data_BNPP.txt','\t',1);
B1=A.data(:,:);

B1(:,2)=B1(:,2)/S01;
S01=1;
B1(:,1)=B1(:,1)/252;
T1=21/252;
B1=B1(B1(:,1)==T1,:);

P1=B1;
for i=1:size(B1,1)
    P1(i,3)=european_bs(S01,B1(i,2),r1,B1(i,3),T1,"call");
end

for i=1:size(B1,1)
      euro1=@(sigma)european_bs(S01,P1(i,2),r1,sigma,T1,"call")-P1(i,3);
      euro2=@(sigma)european_bs(S01,P1(i,2),r1,sigma,T1,"call")-P1(i,3)*0.99998;
      euro3=@(sigma)european_bs(S01,P1(i,2),r1,sigma,T1,"call")-P1(i,3)*1.00001;
      
      C1(i)=fzero(euro1,0.25);
      C2(i)=fzero(euro2,0.5);
      C3(i)=fzero(euro3,0.5);
end

    scatter(P1(:,2),C1(:),'.');
    hold on;
    scatter(P1(:,2),C2(:),'x');
    hold on;
    scatter(P1(:,2),C3(:),'+');
    hold on;
    
    european_bs(S01,0.5025,r1,3.1825e-16,T1,"call")
    
    
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