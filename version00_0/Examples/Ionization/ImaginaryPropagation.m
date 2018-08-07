
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Import the wavefunction spatial axis
x=importdata('out0.txt');
%Import the time axis for the wf plot
t=importdata('out4.txt');
%Import the wavefuntion vs t
A=importdata('out5.txt');
n=800
B=reshape(A,n,length(A)/n);

[T X]=meshgrid(t,x);
figure 
surf(T,X,log10(B+1e-10),'FaceColor','interp',...
 'EdgeColor','none')
title('Evolution of the WF')
xlabel('Time (a.u.)','fontsize',12,'fontweight','b')
ylabel('x (a.u.)','fontsize',12,'fontweight','b')
%daspect([5 5 1])
caxis([-10 0])

axis tight
view(2)

figure
plot(x,B(:,length(A)/n),'r','LineWidth',5)
title('Final WF')
xlabel('x (a.u.)','fontsize',12,'fontweight','b')
ylabel('|phi|^2','fontsize',12,'fontweight','b')


