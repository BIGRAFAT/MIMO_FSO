clc;clear all;close all

tic

load mySISO
load myMIMO
thick=[200 225 250 275 300];
Eb_No=-5:40;

figure
semilogy(Eb_No,SISO_cloud','o-'),xlim([Eb_No(1) Eb_No(end)]),ylim([1e-6 1])
legend('200m','225m','250m','275m','300m')
grid on
title('SISO withCloud noFading')
figure
semilogy(Eb_No,SISO_cloud_fading','o-'),xlim([Eb_No(1) Eb_No(end)]),ylim([1e-6 1])
legend('200m','225m','250m','275m','300m')
grid on
title('SISO withCloud withFading')

figure
semilogy(Eb_No,MIMO_cloud','o-'),xlim([Eb_No(1) Eb_No(end)]),ylim([1e-6 1])
grid on
title('MIMO withCloud noFading')
figure
semilogy(Eb_No,MIMO_cloud_fading','o-'),xlim([Eb_No(1) Eb_No(end)]),ylim([1e-6 1])
grid on
title('MIMO withCloud withFading')

break

for j=1:length(Eb_No)
    figure
    % Fig  4.6
    semilogy(thick,SISO_cloud(:,j),'gx-','LineWidth',2),hold on
    semilogy(thick,MIMO_cloud(:,j),'ro-','LineWidth',2)
    semilogy(thick,SISO_cloud_fading(:,j),'cp-','LineWidth',2)
    semilogy(thick,MIMO_cloud_fading(:,j),'m^-','LineWidth',2),xlim([thick(1) thick(end)]),ylim([1e-6 1])
    legend('SISO Cloud','MIMO Cloud','SISO Cloud Fading','MIMO Cloud Fading')
    grid on
    title(sprintf('at Eb/No = %d dB',Eb_No(j)))
end

toc