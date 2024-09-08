function gril=grille2(xmin,xmax,dx,ymin,ymax,dy)
%
% function gril=grille2(xmin,xmax,dx,ymin,ymax,dy);
%
%  fonction pour creer une grille reguliere 2D
%  allant par pas de dx et dy , a partir de xmin, ymin jusqu'a xmax,ymax
%
x=[xmin:dx:xmax]';
y=[ymin:dy:ymax]';
nx=length(x);
ny=length(y);
gril=[kron(ones(ny,1),x), kron(y,ones(nx,1))];