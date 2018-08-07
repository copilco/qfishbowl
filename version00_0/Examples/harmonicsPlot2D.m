list=importdata('lista')

n1=40
n2=40

figure
hh=1
for i=1:length(list) %4:16%length(list)

    
    a=importdata(list{i});
b=reshape(a,n1,n2);

%subplot(2,2,hh)
surf(log10(b+1e-10),'FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
%daspect([5 5 1])
axis tight
view([90 90])
camlight left
pause(0.2)

end
