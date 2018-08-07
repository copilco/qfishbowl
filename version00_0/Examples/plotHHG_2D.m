


a=importdata('2dout7.txt');
b0=a(:,2);

dt=0.05
dw=2*pi/length(b0)/dt;

c0=fft(b0);
w=-pi/dt:dw:pi/dt-dw;

%subplot(1,2,2)
plot(w/0.057,log10(fftshift(abs(c0))+1e-10) );
xlim([0 80])

w0=0.057
e0=sqrt(2e14/3.5e16)
Up=e0*e0/4/w0/w0;
cutoff=(Up*3.17+0.45)/w0;
line([cutoff cutoff],[-4 4])

title('HHG 2D, Ip=-0.5 a.u. I=2e14 W/cm2 w=0.057 au')
xlabel('Harmonic Oder')
ylabel('HHG Signal')