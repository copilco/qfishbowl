list=importdata('lista')

n1=80
n2=80

figure
hh=0
for i=1:length(list) %4:16%length(list)
hh=hh+1
    
    a=importdata(list{i});
    b=reshape(a,n1,n2);

%subplot(2,2,hh)
surf(log10(b+1e-6),'FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
%daspect([5 5 1])
axis tight
view([90 90])
camlight left
pause(0.2)

end
