
function xx =ICT(Tf,M)
ITERATE_NUM=M;
ofdm_signal=Tf(:);
CR=3;
Signal_Power = abs(ofdm_signal.^2);
Peak_Power = max(Signal_Power);
Mean_Power = mean(Signal_Power);
for nIter=1:ITERATE_NUM
% Clipping
x_tmp = ofdm_signal(Signal_Power>CR*Mean_Power);
x_tmp = sqrt(CR*Mean_Power)*x_tmp./abs(x_tmp);
ofdm_signal(Signal_Power>CR*Mean_Power) = x_tmp;
Signal_Power = abs(ofdm_signal.^2);
% Peak_Power = max(Signal_Power);
% Mean_Power = mean(Signal_Power);
end
xx=reshape(ofdm_signal,size(Tf));
% k=0.5;b=0.01;
% fx=lstcom(k,xx,b);