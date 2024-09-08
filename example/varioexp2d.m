function [gexp]=varioexp2d(x,nbclas,lclas,vdir,vreg)
% Fonction pour calculer plusieurs variogrammes en 2D selon des directions et régularisations spécifiées.
%
% syntaxe function [gexp]=varioexp2d(x,nbclas,lclas,vdir,vreg)
%
%
% input   x : matrice n*3 contenant les coordonnees (x,y) et la variable;
%         nbclas: nombre de classes dans les calculs de vario.
%         lclas: largeur des classes (scalaire)
%                ou matrice(nbclas x 2) donnant les limites des classes (>binf, <=bsup)
%         vdir: vecteur ndirx1.  Directions (en azimut)  
%         vreg: regularisation (degre 0-90) angle de part et d'autre de la direction considérée; vecteur (ndir x 1);
%
% output  gexp: cube contenant les distances moyennes (1 ere colonne)
%                              le nombre de paires    (2e colonne)
%                              le variogramme expérimental(3e colonne)
%                              1 plan par direction demandée et dans l'ordre demandé
%
% auteur D. Marcotte, avril. 2002


if length(lclas)==1
   lclas=[[0:lclas:(nbclas-1)*lclas]',[lclas:lclas:(nbclas)*lclas]'];
end
ncl=length(lclas);

if size(vdir,1)~=length(vreg)
   display('Fournir autant de directions que de régularisations')
   return
end
ndir=size(vdir,1);

n=length(x);
k=0;
var=[];
[n,p]=size(x);
xx=x(:,1);
yy=x(:,2);
v=x(:,3);
gexp=zeros(ncl,3,ndir);
nplot=ceil(sqrt(ndir));
u=poletocart([vdir zeros(ndir,1)]); % u contient les directions pendages sous forme de vecteurs unitaires (cosinus directeurs)
tol=cosd(vreg);

% former les indices de paires


for i=1:n-1,  % boucle des observations
   ind=[ones((n-i),1)*i,[(i+1):n]'];
   yt=ind;                     % vecteur des indices  
   
   % Calcul des distances 
   
   ht=sqrt((xx(yt(:,1))-xx(yt(:,2))).^2+(yy(yt(:,1))-yy(yt(:,2))).^2);
   
   
   % calcul des angles
   dx=xx(yt(:,2))-xx(yt(:,1));
   dy=yy(yt(:,2))-yy(yt(:,1));
  
   uobs=[dx,dy]./(sqrt(dx.^2+dy.^2)*ones(1,2)); % vecteur unitaire correspondant aux paires;
   uobs(:,3)=0;
   
   % trouver les limites angulaires et identifier le cas
   for idir=1:ndir;  % boucle des directions
      
      da=uobs*u(idir,:)';
      id=abs(da)>=tol(idir);
      
      % ne conserver que les paires entrant dans l'angle recherché
      
      h=ht(id);
      y=yt(id,:);
      
      % Calcul des carres des ecarts des variables
      
      var=0.5*((v(y(:,1))-v(y(:,2))).^2);
      
      % Tri sur h et calcul des carres des ecarts correspondants
      
      for ic=1:ncl
         id=h>lclas(ic,1) & h<=lclas(ic,2);
         gexp(ic,2,idir)=gexp(ic,2,idir)+sum(id);
         gexp(ic,1,idir)=gexp(ic,1,idir)+sum(h(id));
         gexp(ic,3,idir)=gexp(ic,3,idir)+sum(var(id));                                         
      end
   end % boucle des directions
   
end % boucle des observations

% calculer le variogramme

for idir=1:ndir;
   id=gexp(:,2,idir)>0;
   gexp(id,1,idir)=gexp(id,1,idir)./gexp(id,2,idir);
   gexp(id,3,idir)=gexp(id,3,idir)./gexp(id,2,idir);
end

ymax=max(max(gexp(:,3,:)))*1.2;
%ymax=max(ymax,1);
xmax=max(max(gexp(:,1,:)))*1.2;
xmax=max(xmax,1);

%for idir=1:ndir   
%   
%   id=gexp(:,2,idir)>0;
%   subplot(nplot,nplot,idir)
%   plotgExp(gexp(id,1,idir),gexp(id,3,idir),gexp(id,2,idir),'+',6)
%   axis([0 xmax, 0, ymax]);
%   title(['Vario. exp. (dir, reg) : (',num2str(vdir(idir,1),3),',',num2str(vreg(idir),3),')']);
%   xlabel('h')
%   ylabel('\gamma(h)')
%end


function x=poletocart(pole)
%
% fonction pour convertir des pôles en coordonnées cartésiennes (système main droite)
% pole: nx2 dir (azimut), pendage du pole angle positif vers le bas (retourne une coordonnée z négative)
deg2rad=pi/180;
pole=pole*deg2rad;

x(:,1)=cos(pole(:,2)).*sin(pole(:,1));
x(:,2)=cos(pole(:,2)).*cos(pole(:,1));
x(:,3)=-sin(pole(:,2));


function s=cosd(a);
% function qui retourne le cosinus de a (en degrés)
a=a/180*pi;
s=cos(a);
