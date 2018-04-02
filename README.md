# Abstract
Free space optical communication (FSO) is one of the sprouting technologies in optical
communication systems domain. It is a cost effective and high bandwidth access technique,
which has been receiving growing attention with recent commercialization successes; thus it
can be employed as an alternative for the conventional radio frequency (RF) links to work out
the current limitations in communication systems. But, FSO from satellite and ground or air to
air consists of cloud as part of the communication channel. Cloud causes pulse broadening in
time domain and attenuation in radiant pulse power. So, the bit error rate (BER) performance
decreases with thickness of cloud. Spatial diversity, which involves the use of multiple laser
transmitters or receivers or both, can be used over FSO links to mitigate cloud-induced
complications with fading channels. Here, we have employed multiple-input multiple-output
(MIMO) technique to improve the BER performance of an FSO system. We have analyzed
BER of single-input single-output (SISO) system with MIMO system undergoing cloud
channels of different thickness using binary phase shift keying (BPSK) signaling technique.
Significant improvement in receiver performance is obtained by simulating the systems under
same conditions.




`thesis report`:
* presents the findings

`siso.m`:
* simulates the transmission & reception of SISO system in Rayleigh Fading & AWGN channel
* produces Figure 4.1 (BER vs Eb/No for SISO with cloud & no fading)
* produces Figure 4.2 (BER vs Eb/No for SISO with cloud & fading)
* saves the obtained values to `mysiso.mat`

`mimo.m`:
* simulates the transmission & reception of MIMO system in Rayleigh Fading & AWGN channel
* produces Figure 4.3 (BER vs Eb/No for MIMO with cloud & no fading)
* produces Figure 4.4 (BER vs Eb/No for MIMO with cloud & fading)
* saves the obtained values to `mymimo.mat`

`powerpenalty.m`:
* produces Figure 4.5 (power penalty vs cloud thickness for various configarations)

`ber.m`:
* produces Figure 4.6 (ber vs cloud thickness for various configarations)

