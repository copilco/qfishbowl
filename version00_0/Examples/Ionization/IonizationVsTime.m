A=importdata('out3.txt');
B=importdata('out6.txt');

e0=max(A(:,2))

hold on
plot(B(:,1),1-B(:,2),'b','LineWidth',5)
scaler=max(1-B(:,2))
plot(A(:,1),A(:,2)*scaler/e0,'r','LineWidth',5)
title('Ionization vs field','fontsize',12,'fontweight','b')
xlabel('Time (a.u.)','fontsize',12,'fontweight','b')
ylabel('Ionization / E (a.u.)','fontsize',12,'fontweight','b')