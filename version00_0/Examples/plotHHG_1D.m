a=importdata('out7.txt');
b0=a(:,2);

dt=0.05
w0=0.057
dw=2*pi/length(b0)/dt;

c0=fft(b0);
w=-pi/dt:dw:pi/dt-dw;

plot(w/w0,log10(fftshift(abs(c0))+1e-10) );
xlim([0 80])

e0=sqrt(2e14/3.5e16)
Up=e0*e0/4/w0/w0;
cutoff=(Up*3.17+0.58)/w0;
line([cutoff cutoff],[-2.5 2.5])
title('HHG 1D, Ip=-0.5 a.u. I=2e14 W/cm2 w=0.057 au')
xlabel('Harmonic Oder')
ylabel('HHG Signal')