function gril=grille3(xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz);
%
% function gril=grille3(xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz);
%
%  fonction pour creer une grille reguliere 3D
%  allant par pas de dx dy et dz, a partir de xmin ymin et zmin
%  jusqu'a xmax,ymax,zmax
%
x=[xmin:dx:xmax]';
y=[ymin:dy:ymax]';
z=[zmin:dz:zmax]';
nx=length(x);
ny=length(y);
nz=length(z);
gril=[kron(ones(ny*nz,1),x), kron(kron(ones(nz,1),y),ones(nx,1)), kron(z,ones(nx*ny,1))];
end