function [datasim,S,u,G,GU,U,z]=fftma(model,c,seed,nbsimul,nx,dx,ny,dy,nz,dz)
% fftma: fonction pour simuler par FFT-MA (voir Le Ravalec, Math. Geol. 2000, #32, 701-723)
% syntaxe: z=fftma(model,c,seed,nbsimul,nx,dx,ny,dy,nz,dz)
% Simulation par FFT-MA
% si 1D spécifier uniquement nx,dx,
% si 2D spécifier nx,dx,ny,dy
% si 3d, spécifier nx,dx,ny,dy,nz,dz
% input: Nx,dx,Ny,dy,Nz,dz: nombre et pas de la grille en x,y,z 
% 
% model: modele de variogramme (comme cokri)
% c: plateaux (comme cokri)
% seed, germe aléatoire
% nbsimul: nombre de réalisations
%
% ATTENTION: le champ a simuler est (nx+a)*(ny+a)*(nz+a) ou a est la portée
% effective maximale des différentes composantes du modele.
% Éviter de spécifier un modele avec une tres grande portée. 

% auteurs D. Marcotte et E. Henry oct 2003

if nargin==6 
    cas=1;
    ny=1;nz=1;
    dy=1;dz=1;
elseif nargin==8
    cas=2;
    nz=1;dz=1;
elseif nargin==10
    cas=3;
else
    errordlg('Vérifiez l''input, nombre de parametres incorrect')
end
rng('default')
rng(seed);
%randn('state',seed);

[Nx,Ny,Nz]=nombremin(model,c,nx,dx,ny,dy,nz,dz,cas);

Nx2=floor(Nx/2);
Ny2=floor(Ny/2);
Nz2=floor(Nz/2);

% Echantillonnage de la covariance C
switch cas
    case 1
        C=ones(Nx,1)*nan;
        x0=[-(Nx2+1)*dx:dx:Nx2*dx]';
        cc=covar(x0,0,model,c);
        C(1:Nx+1)=cc;
        
    case 2
        C=ones(Nx,Ny)*nan;
        x0=grille2(-(Nx2+1)*dx,Nx2*dx,dx,-(Ny2+1)*dy,Ny2*dy,dy);
        cc=covar(x0,[0 0],model,c);
        cc=reshape(cc,Nx+1,Ny+1);
        C(1:Nx+1,1:Ny+1)=cc;
        
    case 3
        
        C=ones(Nx,Ny,Nz)*nan;
        x0=grille3(-(Nx2+1)*dx,Nx2*dx,dx,-(Ny2+1)*dy,Ny2*dy,dy,-(Nz2+1)*dz,Nz2*dz,dz);
        cc=covar(x0,[0 0 0],model,c);
        cc=reshape(cc,Nx+1,Ny+1,Nz+1);
        C(1:Nx+1,1:Ny+1,1:Nz+1)=cc;
        
end

C=fftshift(C);

% Génération des valeurs u

S=fftn(C);

% Calcul de G

G=sqrt(S);


% définition des coordonnées
switch cas
    case 1
        x0=[0:dx:(nx-1)*dx]';
    case 2
        x0=grille2(0,(nx-1)*dx,dx,0,(ny-1)*dy,dy);
    case 3
        x0=grille3(0,(nx-1)*dx,dx,0,(ny-1)*dy,dy,0,(nz-1)*dz,dz);
end

% simulation des différentes réalisations
datasim = [];
for i=1:nbsimul
    
    u=randn(size(C));    
    U=fftn(u);    
    GU=G.*U;    
    % Transformation de Fourier inverse donnant g*u et z    
    z=real(ifftn(GU));    
    z=z(1:nx,1:ny,1:nz);
    z=reshape(z,nx*ny*nz,1);    
    datasim = [datasim z]; 
end




function [nxsim,nysim,nzsim]=nombremin(model,c,nx,dx,ny,dy,nz,dz,cas)
% fonction pour determiner le nombre de points minimum a simuler ce nombre
% prend en compte la portée du a la périodicite de la FFT
fact=1;
id=model(:,1)==1;
model=model(~id,:);
ctot=sum(c(~id));
[nm,p]=size(model);

if p==2
    if cas == 2
        model=[ model(:,:) model(:,2) zeros(nm,1) ];
        p=4;
    end
end

a=model(:,2:min(p,p-1+cas));

for i=1:nm
    if model(i)==2
        a(i,1:cas)=a(i,1:cas)*3;
    elseif model(i)==3;
        a(i,1:cas)=a(i,1:cas)*sqrt(3);
    elseif model(i)==8;
        a(i,1:cas)=a(i,1:cas)*20;
    elseif model(i)==9;
        a(i,1:cas)=a(i,1:cas)*2.6;
    elseif model(i)==19;
        a(i,1:cas)=a(i,1:cas)*3;
    end
end

if cas==2
    a(max(a(:,1:2),[],2)~=a(:,1),3)=a(max(a(:,1:2),[],2)~=a(:,1),3)+90;
    a(max(a(:,1:2),[],2)~=a(:,1),[1,2])=a(max(a(:,1:2),[],2)~=a(:,1),[2,1]);
    
    ag=a(:,1);
    ap=a(:,2);
    thetax=360-a(:,3);
    thetay=450-a(:,3);
    
    ax=(ap.*ag)./sqrt( (ap.^2).*(cos(thetax*pi/180).^2) + (ag.^2).*(sin(thetax*pi/180).^2) );
    ay=(ap.*ag)./sqrt( (ap.^2).*(cos(thetay*pi/180).^2) + (ag.^2).*(sin(thetay*pi/180).^2) );
    
    
    clear a
    a(:,1:cas)=max([ax,ay]);
end
amaxx=a(1,1);
nxsim=ceil(fact*amaxx/dx)+1+nx;
nysim=1;
nzsim=1;
if cas==2  % on est en 2D
    amaxy=a(1,2);
    nysim=ceil(fact*amaxy/dy)+1+ny;
    if ny==1 % dans ce cas on est en realité en 1D
        nysim=1;
    end
    nzsim=1;
elseif cas==3
    amaxy=a(1,2);
    nysim=ceil(fact*amaxy/dy)+1+ny;
    if ny==1
        nysim=1;
    end
    amaxz=a(1,3);
    nzsim=ceil(fact*amaxz/dz)+1+nz;
    if nz==1
        nzsim=1;
    end
end

if mod(nxsim,2)==0
    nxsim=nxsim+1;
end
if mod(nysim,2)==0
    nysim=nysim+1;
end
if mod(nzsim,2)==0
    nzsim=nzsim+1;
end
