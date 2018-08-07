%% READER BOUND STATES FOR HYDROGEN ATOM

clear all

A0=importdata('out0.txt');
S0=importdata('state0.txt');
S1=importdata('state1.txt');
S2=importdata('state2.txt');

%% Soft-Core Hydrogen potential well 
V=@(x,a) -1.0./sqrt(a+x.^2);

%% Loading axis
x   = A0(:,1);
n   = length(x);



%% Potential well representation %%

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)...
    scrsz(3)*0.85 scrsz(4)*0.8],'Color','w')

plot(x,V(x,2),'k','linewidth',2)
hold on
plot(x,(S0(2:n+1,1)+1e-18).^2+S0(1,1),'b')
%plot(x,log10(A2(1+nx*(ktime-1):ktime*nx)+1e-12),'b')
%plot(x,log10(A3(1+nx*(ktime-1):ktime*nx)+1e-12),'b')

axis tight
grid off


ylabel('  |\psi|^2 ','fontsize',16,'fontweight','b')
xlabel(' x (a.u.) ','fontsize',16,'fontweight','b')

title(['Energy of the ground state; E0: ',num2str(S0(1,1)),' a.u. ']...
    ,'fontsize',16,'fontweight','b')


xlim([-30 30])
set(gca,'fontsize',16,'fontweight','b')
box off


%% End visualization %%  