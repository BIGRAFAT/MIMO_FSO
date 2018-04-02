% MIMO
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

fc=1e6;
Tc=1/fc;
T=1*Tc; %%%%%%%%%%%%%%%%
% Tc=T/10;
Ts=T/b;
t=0:Ts:T-Ts;
sin_fn=sin(2*pi*fc*t);
N=n*b;
f=(-N/2+1:N/2)*fc/n;
Dr=.25;


D=[ 1 1; 1 0; 0 1; 0 0];
h_mimo=[.8 .2;.5 -.5];


flag=2*sin_fn(1)-1;
flag=[1 -flag;1 flag];
mod_fn=kron(flag,sin_fn);
demod_fn=mod_fn;

lpf_demod=zeros(2,b);
lpf_demod(:,1:round(2*T/Tc))=1; %%%%%%%%%%%%%%%    "2" is used in LPF
lpf_demod=[lpf_demod fliplr(lpf_demod)];

D_mod=zeros(2,2*b,4);
D_demod=D_mod;
for j=1:4
    d=D(j,:);
    d_map=zeros(2,2);
    d_map(:,1:2:end)=reshape(d,2,[]); % [x1 x2  ...]
    d_map(:,2:2:end)=flipud(d_map(:,1:2:end)); % [-x2* x1* ....]
    d_map=2*d_map-1;
    d_mod=kron(d_map,ones(1,b)).*mod_fn; % TX_BPSK
    d_demod=d_mod.*demod_fn;
    D_mod(:,:,j)=d_mod;
    D_demod(:,:,j)=d_demod;
end

tx=zeros(2,n*b);
for j=1:2:n
    temp=in(j:j+1);
    q=4-( 2*temp(:,1)+temp(:,2) );
    tx( : , (j-1)*b+1 : (j+1)*b )=D_mod(:,:,q);
end

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

h_rayleigh = ( randn(2,n*b) + 1i*randn(2,n*b) ) / sqrt(2) ; % Rayleigh channel
% h_rician = ( (k*cos(t)+randn(1,n)) + 1i*(k*sin(t)+randn(1,n)) ) / sqrt(2) ; % Rician channel

ber=zeros(5,length(Eb_No),2);
er=zeros(1,4);
for c=1:5
    c
    H=H_all(c,:);
    h=h_all(c,:);
    
    tx_conv=[];
    tx_conv(1,:)=conv(tx(1,:),h);
    tx_conv(2,:)=conv(tx(2,:),h);
    tx_conv=tx_conv(:,1:length(tx));
    
    tx_mix_all(1:2,:)=h_mimo*tx_conv;
    
    tx_conv=tx_conv.*h_rayleigh;
    tx_mix_all(3:4,:)=h_mimo*tx_conv;

    
    for k=1:2
        tx_mix=tx_mix_all(2*k-1:2*k,:);
        for j=1:length(Eb_No)
            s=Eb_No(j);
            ber_sum=0;
            for p=1:5
                noise=( randn(1,n*b) + 1i*randn(1,n*b) ) / sqrt(2); % white gaussian noise, 0dB variance
                rx_mix(1,:)=tx_mix(1,:)+noise*10^(-s/10)*sqrt(1/2)*(1/2);
                noise=( randn(1,n*b) + 1i*randn(1,n*b) ) / sqrt(2); % white gaussian noise, 0dB variance
                rx_mix(2,:)=tx_mix(2,:)+noise*10^(-s/10)*(1/2)*(1/2);

    %             rx_mix=rx_mix./h_rayleigh;
                rx=h_mimo\rx_mix;
                if k==2
                    rx=rx./h_rayleigh;
                end
                rx(1,:)=ifft( fft(rx(1,:))./H) ;   % EQUALIZER filter = Hc(f) = 1/H(f)
                rx(2,:)=ifft( fft(rx(2,:))./H) ;   % EQUALIZER filter = Hc(f) = 1/H(f)

                out=zeros(1,n);
                for v=1:2:n
                    temp=rx( : , (v-1)*b+1 : (v+1)*b );
                    temp_demod=temp.*demod_fn;
                    fft_temp_demod=fft(temp_demod');
                    lpf_temp_demod=real(ifft(fft_temp_demod.*lpf_demod'));    %  LPF to smoothen
                    for w=1:4
                        d_demod=D_demod(:,:,w);
                        er(w)=sum(sum( (d_demod-lpf_temp_demod').^2 ));
                    end
                    flag_out=find( er==min(er) , 1 );
                    out(v:v+1)=D(flag_out,:);
                end
                ber_sum=ber_sum+sum(out~=in);
            end
            ber(c,j,k)=ber_sum/p/n;
        end
        toc
    end
    
end
MIMO_cloud=ber(:,:,1);
MIMO_cloud_fading=ber(:,:,2);

figure
% Fig. 4.3
semilogy(Eb_No,MIMO_cloud','o-'),xlim([Eb_No(1) Eb_No(end)]),ylim([1e-6 1])
grid on
title('MIMO withCloud noFading')

figure
% Fig. 4.4
semilogy(Eb_No,MIMO_cloud_fading','o-'),xlim([Eb_No(1) Eb_No(end)]),ylim([1e-6 1])
grid on
title('MIMO withCloud withFading')


save MIMO MIMO_*

