% SISO
% withCloud noFading
% withCloud withFading

clc;clear all;
% close all

tic
color='bgrcmyk';
thick=[200 225 250 275 300];
k1=[120.1   34.1    12.4    5.1     2.4];
k2=[1.9     1.9     1.1     0.8     1.7]*1e7;
k3=[1.55    1.6     0.66    0.28    0.19];
k4=[3       3       2.4     1.8     1.6]*1e6;

n=1000; %%%%%%%%%%%%
in=randsrc(1,n,[0 1]);
b=10; %%%%%%%%%%%%%%%
Eb_No=-5:40;

fc=10e6;
Tc=1/fc;
T=1*Tc; %%%%%%%%%%%%%%%%
% Tc=T/10;
Ts=T/b;
t=0:Ts:T-Ts;
sin_fn=sin(2*pi*fc*t);
N=n*b
f=(-N/2+1:N/2)*fc/n;
Dr=.20;
demod_fn=kron(ones(1,n),sin_fn);
tx=kron(2*in-1,sin_fn); % TX_BPSK
TX=fft(tx);

h_rayleigh = ( randn(1,n*b) + 1i*randn(1,n*b) ) / sqrt(2) ; % Rayleigh channel
% h_rician = ( (k*cos(t)+randn(1,n)) + 1i*(k*sin(t)+randn(1,n)) ) / sqrt(2) ; % Rician channel

h_all=zeros(5,n*b);
H_all=h_all;
for c=1:5
    H=k1(c)./(k2(c)+1i*2*pi*f).^2 + k3(c)./(k4(c)+1i*2*pi*f).^2;
    H=H*pi*Dr^2/4*fc^2;
    H=ifftshift(H);
    h=ifft(H);
    h_all(c,:)=h;
    H_all(c,:)=H;
    
end

ber=zeros(5,length(Eb_No),2);
for c=1:5
c
    H=H_all(c,:);
    h=h_all(c,:);

%     tx=tx.*h_rayleigh;
    y_conv=conv(tx,h);
    y_conv=y_conv(1:length(tx));
    y_all(1,:)=y_conv;
    y_all(2,:)=y_conv.*h_rayleigh;
    
    for k=1:2
        y_now=y_all(k,:);
        for j=1:length(Eb_No)
            s=Eb_No(j);
            ber_sum=0;
            for p=1:5
                noise=1/sqrt(2)*( randn(1,n*b) + 1i*randn(1,n*b) ); % white gaussian noise, 0dB variance
                rx=y_now+noise*10^(-s/10)*sqrt(1/2);
                if k==2
                    rx=rx./h_rayleigh;
                end
                rx=ifft( fft(rx) .* (1./H) );   % EQUALIZER filter = Hc(f) = 1/H(f)
    %             rx=rx./h_rayleigh;
                rx_demod=rx.*demod_fn;
                rx_demod=reshape(rx_demod',b,[])';
                rx_demod_int=sum(rx_demod,2);
                rx_demod_int=reshape(rx_demod_int,[],1)'/1;
                out=rx_demod_int>0;
                ber_sum=ber_sum+sum(in~=out);
            end
            ber(c,j,k)=ber_sum/p/n;
    %         [ s 10*log10( sum(abs(tx)) / (sum(abs(noise))*10^(-s/10)*sqrt(1/2)) )]
        end
        toc
    end
end

SISO_cloud=ber(:,:,1);
SISO_cloud_fading=ber(:,:,2);


figure
% Fig. 4.1
semilogy(Eb_No,SISO_cloud','o-'),xlim([Eb_No(1) Eb_No(end)]),ylim([1e-6 1])
grid on
title('SISO withCloud noFading')

figure
% Fig. 4.2
semilogy(Eb_No,SISO_cloud_fading','o-'),xlim([Eb_No(1) Eb_No(end)]),ylim([1e-6 1])
grid on
title('SISO withCloud withFading')

save mySISO SISO_*

