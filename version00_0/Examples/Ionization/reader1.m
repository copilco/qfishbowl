
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
surf(T,X,log10(B+1e-6),'FaceColor','interp',...
 'EdgeColor','none')
title('Evolution of the WF')
xlabel('Time (a.u.)','fontsize',12,'fontweight','b')
ylabel('x (a.u.)','fontsize',12,'fontweight','b')
%daspect([5 5 1])
caxis([-6 -2])

axis tight
view(2)



C=importdata('out3.txt');
v1=ones(1,length(C));
e0=max(C(:,2))
hold on 
plot3(C(:,1),C(:,2)*50/e0,v1*-2,'r','LineWidth',4);