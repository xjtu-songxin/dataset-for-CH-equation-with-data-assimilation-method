close all;clc;clear 
rng(1);
Nx =100;
Ny = Nx;
xa = 0; xb = 1;
ya = xa ; yb = xb;
h=(xb-xa)/Nx;
x = linspace(xa+0.5*h,xb-0.5*h,Nx);
y = linspace(ya+0.5*h,yb-0.5*h,Ny);
dt = 1e-7;    
epsilon = 0.01;
maxiter = 1e7;
ns=1e6;
nrecord = 1000;
alpha= 0;
beta = 1;
lambda = 10;
theshold = 0;
file_name = 'rate_'+string(theshold)+'_2d.eps';

xaxis = [];
yaxis = []; 

[yy,xx] = meshgrid(y,x);
kx = 2*pi/(xb-xa)*[0:Nx/2 -Nx/2+1:-1];
ky = 2*pi/(yb-ya)*[0:Ny/2 -Ny/2+1:-1];
k2x = kx.^2;   k2y = ky.^2;
k4x = kx.^4;   k4y = ky.^4;
[kxx,kyy]=meshgrid(k2y,k2x);
[kxxxx,kyyyy]=meshgrid(k4x,k4y);
wei=(kxx+kyy);
wei1=(kxxxx + kyyyy + 2.*kxx.*kyy);
 
ave=0.00;
phi=2*rand(Nx,Ny)-1;


figure(1)
contourf(phi)
axis image
colorbar
caxis([-1, 1])

rng(2);
mask = rand(Nx,Ny);
mask(mask>theshold) = 1;
mask(mask<=theshold) = 0;

%% load measurement value
psi=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        if i<Nx/2
            psi(i,j) = 1;
        else
            psi(i,j) = -1;
        end
    end
end
size_psi = size(psi);
interval = maxiter+5;
 
ophi=phi
weight_start = sum(sum(phi))*h*h;
xaxies = []; 
massaxies = []; 
energy = []; 
for it = 1:maxiter
     it
     
     hphi=2.0*phi-ophi;
     
     phif = phi.^3 - 1.0*phi; phi0f = ophi.^3 - 1.0*ophi;

     A = 1 + 0.5*(epsilon^2).*wei1.*dt + 0.5*(lambda.*wei).*dt;
     B = fftn(phi) - 0.5*(epsilon^2).*wei1.*dt.*fftn(phi) - 0.5*beta*dt.*wei.*(fftn(phif) - fftn(phi0f)) ...
         - dt*beta*wei.*fftn(phif)- 0.5*dt*lambda.*wei.*(fftn(-2*phi + ophi)) ...
         - 0.5*wei.*dt*alpha .* fftn(mask.*(hphi + phi-2*psi));
     
     
     phinew = B./A;
     
     ophi = phi;
     
     phi = real(ifftn(phinew));
     
     if (mod(it,nrecord) == 0)
        c_T = it*dt;
        file_name = './reference_solution/'+string(c_T)+'.mat';
        save(file_name, "phi")
     end     
 end
weight_end = sum(sum(phi))*h*h;