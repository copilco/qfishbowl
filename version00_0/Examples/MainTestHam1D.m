
n1=400
n2=200 

a0=importdata('out5.txt');
    b0=reshape(a0,n1,n2);
    surfc(log10(b0+1e-4),'FaceColor','interp',...
        'EdgeColor','none')
        
    axis ([ 0 n2 0 n1 ]);
    grid off
    axis on
    pause(0.2)
end