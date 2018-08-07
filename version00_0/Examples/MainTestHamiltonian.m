lista=importdata('lista')

n1=80
n2=80 
for h=1:length(lista)
    a0=importdata(lista{h});
    b0=reshape(a0,n1,n2);
    surfc(b0,'FaceColor','interp',...
        'EdgeColor','none')
        
    axis([ 0 n1 0 n1 ]);
    grid off
    axis on
    pause(0.2)
end