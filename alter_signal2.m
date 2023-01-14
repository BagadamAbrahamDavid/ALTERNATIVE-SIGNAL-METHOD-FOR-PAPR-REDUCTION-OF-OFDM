clc;
clear all;
close all;
Nsym=2048;
Nsub=64;
L=4;
M=4;
rand('state',1);
r=randi([0 M-1],Nsym,Nsub);
dm=qammod(r,M);
%-----------------------------------------
for i=1:Nsym 
     time_domain_signal=(ifft([dm(i,1:Nsub/2) zeros(1,(L-1)*Nsub) dm(i,Nsub/2+1:Nsub)]));
     meano=mean(abs(time_domain_signal).^2);
     peako=max(abs(time_domain_signal).^2);
     papro(i)=10*log10(peako/meano);    
end
[y,range]=ccdf(papro,0.5);
semilogy(range,y,'-x','Linewidth',2);hold on;
%=================================================
LL=L*Nsub-1;
fc= zeros(1,LL);
Pc(1)=1;
Pc(2)=0.97195983;
Pc(3)=1/sqrt(2);
Pc(4)=(1-Pc(2)^2)^0.5;
for n=1:LL
    PP=0;
    for m=1:3
        PP=PP+2*(-1)^m*Pc(m+1)*cos(2*pi*m*n/L/Nsym);
    end
    fc(n)=Pc(1)+PP;
end
fc=fc/(sum(fc.^2))^0.5;
%B=circshift(fc,Nsub/2)
%=============================================
for i=1:Nsym 
     time_domain_signal=ifft([dm(i,1:Nsub/2) zeros(1,(L-1)*Nsub) dm(i,Nsub/2+1:Nsub)]);
     t_real=real(time_domain_signal);
     t_imag=imag(time_domain_signal);   
     sig=t_real+j*conv(t_imag,fc,'same');             
     meann=mean(abs(sig).^2);
     peakn=max(abs(sig).^2);
     paprn(i)=10*log10(peakn/meann);    
end
[y1,range1]=ccdf(paprn,0.5);
semilogy(range1,y1,'k-','Linewidth',2);hold on;
ylim([10^-3 10^0]);
%============================================
%======== Applying SLM approach =============
randn('state',1);
C=8;p=[1 -1 j -j]; 
B=randsrc(C,Nsub,p); 
for ii=1:Nsym
P{1}=[dm(ii,1:Nsub/L) zeros(1,3*(Nsub/L))];
V=2;
  for jj=V:L
      if V~=0
P{jj}=[zeros(1,(jj-1)*(Nsub/L)) dm(ii,((jj-1)*(Nsub/L)+1):(jj)*(Nsub/L)) zeros(1,V*(Nsub/L))];
      else
P{jj}=[zeros(1,(jj-1)*(Nsub/L)) dm(ii,((jj-1)*(Nsub/L)+1):(jj)*(Nsub/L))];
      end          
V=V-1;     
  end
%------Transform Pi to Time Domain -------------
 for kk=1:L
     PP=P{kk};
     tmp=(ifft([PP(1:(Nsub/2)) zeros(1,(L-1)*Nsub) PP((Nsub/2+1):Nsub)]));   
     t_real=real(tmp);
     t_imag=imag(tmp);
     pt{kk}=t_real+j*conv(t_imag,fc,'same');           
 end  
%----------------------------------------------
 for kk=1:size(B,1)
     bb=B(kk,:);
     r=pt{1}*bb(1)+pt{2}*bb(2)+pt{3}*bb(3)+pt{4}*bb(4);
     mm(kk)=max(abs(r).^2);
     rst{kk}=r;                       
 end
 %--- conventional system -----------------
 ov=min(mm);opt=find(ov==mm);
 xcopt=rst{opt}; 
 meank=mean(abs(xcopt).^2);
 peak=max(abs(xcopt).^2);
 papr(ii)=10*log10(peak/meank); 
 end 
 [cy,cx]=ccdf(papr,0.5);
 semilogy(cx,cy,'m','Linewidth',2) ;hold on; 
 grid on;
 xlabel('PAPR(dB)');
 ylabel('-CCDF');
%============================================================
% randn('state',0);
% C=8;p=[1 -1 j -j]; 
% B=randsrc(C,Nsub,p); 
% clear pt
% for ii=1:Nsym
% P{1}=[dm(ii,1:Nsub/L) zeros(1,3*(Nsub/L))];
% V=2;
%   for jj=V:L
%       if V~=0
% P{jj}=[zeros(1,(jj-1)*(Nsub/L)) dm(ii,((jj-1)*(Nsub/L)+1):(jj)*(Nsub/L)) zeros(1,V*(Nsub/L))];
%       else
% P{jj}=[zeros(1,(jj-1)*(Nsub/L)) dm(ii,((jj-1)*(Nsub/L)+1):(jj)*(Nsub/L))];
%       end          
% V=V-1;     
%   end
% %------Transform Pi to Time Domain -------------
%  for kk=1:L
%      PP=P{kk};
%      tmp=(ifft([PP(1:(Nsub/2)) zeros(1,(L-1)*Nsub) PP((Nsub/2+1):Nsub)]));   
%      t_real=real(tmp);
%      t_imag=imag(tmp);
%      pt{kk}=t_real+j*conv(t_imag,fc,'same');           
%  end  
% %----------------------------------------------
%  for kk=1:C
%      bb=B(kk,:);
%      R=0;it=1;
%      for pp=1:4
%      r=0;     
%      for tt=1:4         
%      r=r+pt{pp}*bb(it);
%      it=it+1;
%      end
%      R=R+r;
%      end
%      mm=max(abs(R).^2);
%      rst{kk}=R;                       
%  end
%  %--- conventional system -----------------
%  ov=min(mm);opt=find(ov==mm);
%  xcopt=rst{opt}; 
%  meank=mean(abs(xcopt).^2);
%  peak=max(abs(xcopt).^2);
%  paprk(ii)=abs(papro(ii)-10*log10(peak/meank));
%  
%  end 
%  [cyk,cxk]=ccdf(paprk,0.5);
%  semilogy(cxk,cyk,'r','Linewidth',2) ;hold on; 
%  grid on;
%  xlabel('PAPR(dB)');
%  ylabel('-CCDF');
%  legend('Original','OSLM','AS-I','AS-J');