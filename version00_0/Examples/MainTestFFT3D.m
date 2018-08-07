a0=importdata('fileinX.txt');
a1=importdata('fileinQ.txt');
a2=importdata('fileoutX.txt');

b0=reshape(a0,40,40,40);
b1=reshape(a1,40,40,40);
bb1=fftshift(b1);
b2=reshape(a2,40,40,40);

const=1e-3;
const2=1e-3;

subplot (1,3,1)
p = patch(isosurface(b0,const));
isonormals(b0,p)
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1])
view(3); axis tight
camlight 
lighting gouraud

subplot (1,3,2)
p1 = patch(isosurface(bb1,const2));
isonormals(bb1,p1)
set(p1,'FaceColor','green','EdgeColor','none');
daspect([1 1 1])
view(3); axis tight
camlight 
lighting gouraud

subplot (1,3,3)
p2 = patch(isosurface(b2,const));
isonormals(b2,p2)
set(p2,'FaceColor','red','EdgeColor','none');
daspect([1 1 1])
view(3); axis tight
camlight 
lighting gouraud