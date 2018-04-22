clear;
S0=167.65;
r=0.06;

A = importdata('AAPL.txt','\t',1);
B=A.data(:,:);


for idx=unique(B(:,1))'
    
C=zeros(size(B(B(:,1)==idx,3),1),1);
C(:,1)=B(B(:,1)==idx,3);
D=zeros(size(B(B(:,1)==idx,3),1),1);
D(:,1)=B(B(:,1)==idx,3);
A=zeros(size(B(B(:,1)==idx,3),1),1);
A(:,1)=B(B(:,1)==idx,3);

while 1
    bk=0;
for i=2:size(B(B(:,1)==idx,3),1)
    if C(i,1)>C(i-1,1)
        C(i-1,1)=C(i,1);
        bk=1;
    end
    if D(i,1)>D(i-1,1)
        D(i,1)=D(i-1,1);
        bk=1;
    end
A(i,1)=(C(i,1)+D(i,1))/2;
end
if bk==0
    break
end
end
B(B(:,1)==idx,3)=A;
end






%B=B(B(:,2)>Minim & B(:,2)<Maxim   & B(:,1)>Exc & B(:,1)<Excm,:);


%{
for i=1:size(B,1)
euro=@(sigma)european_bs(S0,B(i,2),r,sigma,B(i,1)./365,'call')-B(i,3);
B(i,3)=fzero(euro,0.5);
end
%}
%B=B(B(:,1)==21,:);

xlabel('time(days)');
ylabel('strike');
zlabel('price');
scatter3(B(:,1),B(:,2),B(:,3),'.');

%scatter(B(:,2),B(:,3),'filled');
%hold on;
%scatter(A(:,1),C(:,1),'filled');
%hold on;
%scatter(A(:,1),D(:,1),'filled');
%hold on;
%scatter(A(:,1),A(:,2),'filled');


function euro=european_bs(S0,K,r,sigma0,T,putcall)
d1 = (log(S0/K) + (r + 0.5*sigma0^2)*T)/(sigma0*sqrt(T));
d2 = d1 - sigma0*sqrt(T);
N1 = normcdf(d1);
N2 = normcdf(d2);
if putcall=='call'
    euro = S0*N1 - K*exp(-r*T)*N2;
elseif putcall=='put'
    euro = S0*N1 - K*exp(-r*T)*N2 + K*exp(-r*T) - S0;
end
end