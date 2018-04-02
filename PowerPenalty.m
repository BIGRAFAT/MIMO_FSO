
clc;clear all
% close all

color='bgrcmyk';
thick=[200 225 250 275 300];


power_SISO_cloud=[8.1 13.4 15.85 19.15 24.5]
power_SISO_cloud_fading=[16.8 22.0 24.5 27.7 33.2]

power_MIMO_cloud=[5.0 10.2 12.6 15.8 21.25]
power_MIMO_cloud_fading=[13.6 18.9 21.05 24.4 29.95]

penalty_SISO_fading=power_SISO_cloud_fading-power_SISO_cloud
penalty_MIMO_fading=power_MIMO_cloud_fading-power_MIMO_cloud

penalty_SISO_MIMO=power_SISO_cloud-power_MIMO_cloud
penalty_SISO_MIMO_fading=power_SISO_cloud_fading-power_MIMO_cloud_fading

for i=1:5
    a(i)=sum(penalty_SISO_fading(1:i));
    b(i)=sum(penalty_MIMO_fading(1:i));
    aa(i)=sum(penalty_SISO_MIMO(1:i));
    bb(i)=sum(penalty_SISO_MIMO_fading(1:i));
end

figure(1)
subplot(221),plot(thick,a,'o-'),title('SISO Cloud ~ SISO Cloud Fading'),xlabel('Cloud Thickness  (m)'),ylabel('Power Penalty  (dB)')
subplot(222),plot(thick,b,'o-'),title('MIMO Cloud ~ MIMO Cloud Fading'),xlabel('Cloud Thickness  (m)'),ylabel('Power Penalty  (dB)')
subplot(223),plot(thick,aa,'o-'),title('SISO Cloud ~ MIMO Cloud'),xlabel('Cloud Thickness  (m)'),ylabel('Power Penalty  (dB)')
subplot(224),plot(thick,bb,'o-'),title('SISO Cloud Fading ~ MIMO Cloud Fading'),xlabel('Cloud Thickness  (m)'),ylabel('Power Penalty  (dB)')

figure(2)
% Fig  4.5
subplot(221),plot(thick,a,'o-'),title('SISO Cloud ~ SISO Cloud Fading'),ylim([5 10])
subplot(222),plot(thick,b,'o-'),title('MIMO Cloud ~ MIMO Cloud Fading'),ylim([5 10])
subplot(223),plot(thick,aa,'o-'),title('SISO Cloud ~ MIMO Cloud'),ylim([0 5])
subplot(224),plot(thick,bb,'o-'),title('SISO Cloud Fading ~ MIMO Cloud Fading'),ylim([0 5])


