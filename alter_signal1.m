clc;
clear all;
close all;
Nsym=256;    
Nsub=64;
M=4;
K=4;
t=1;
T=1;
z=exp(t);
%===================================
yy=zeros(2.^8,64);
for number=1:2^8   
    S=0;
    for k=1:8
        temp=mod((number-1-S)/2^(k-1),2);
        yy(number,(k-1)*8+1:k*8)=-2*temp+1;     
        S=S+temp*2^(k-1);
    end
end
%======================================
D=randi([0 M-1],Nsub,Nsym);
DM=qammod(D,M);
%================================
L=K*Nsym-1;
filter= zeros(1,L-1);
P(1)=1;
P(2)=0.97195983;
P(3)=1/sqrt(2);
P(4)=(1-P(2)^2)^0.5;
for n=1:L
    PP=0;
    for m=1:3
        PP=PP+2*(-1)^m*P(m+1)*cos(2*pi*m*n/K/Nsym);
    end
    filter(n)=P(1)+PP;
end
filter=filter/(sum(filter.^2))^0.5;
filter=filter*2;
%====================================
DU=(upsample(DM',K))';
D_real=[real(DU),zeros(Nsub,Nsym/2)];
D_imag=imag(DU);
D_imag_shift=[zeros(Nsub,Nsym/2),D_imag(:,1:length(D_imag))];
D_real_filtered=zeros(Nsub,length(D_real)+length(filter)-1);
D_imag_filtered=zeros(Nsub,length(D_imag_shift)+length(filter)-1);
for m=1:Nsub
        signal_real_filtered(m,:)=conv(D_real(m,:),filter);
        signal_imag_filtered(m,:)=conv(D_imag_shift(m,:),filter);
        signal_filtered=signal_real_filtered(m,:)+1i*signal_imag_filtered(m,:);
        Tx=ifft(signal_filtered);%.*exp(1i*(0:Nsub-1).'*((1:length(signal_filtered))*(2*pi/T)/Nsym+pi/2));
        meano=mean(abs(Tx).^2);
peako=max(abs(Tx).^2);
papro(m)=10*log10(peako/meano);  
end

[y4,range4]=ccdf(papro,0.2);
semilogy(range4,y4,'c-','Linewidth',2);hold on;
%======================================