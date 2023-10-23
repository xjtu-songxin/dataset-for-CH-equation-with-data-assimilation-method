close all;clc;clear 
rng(1);
% Nx = 16;
for Nx =[16,32,64,128,256,512,1024]
Ny = Nx;
xa = 0; xb = 1;
ya = xa ; yb = xb;
h=(xb-xa)/Nx;
x = linspace(xa+0.5*h,xb-0.5*h,Nx);
y = linspace(ya+0.5*h,yb-0.5*h,Ny);
dt = 0.1*h;    
epsilon = 0.01;
maxiter = 0.25*0.1/dt;
ns=100;
nrecord = [100000];
file_name = 'with_data_assimilation';
alpha= 5;
beta = 1;
lambda = 50;

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

phi = zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        if i>Nx/2
            phi(i,j) = 1;
        else
            phi(i,j) = -1;
        end
    end
end


rng(2);
mask = rand(Nx,Ny);
theshold = 0;
mask(mask>theshold) = 1;
mask(mask<=theshold) = 0;

%% load measurement value
psi=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        psi(i,j)=tanh((0.45^2-(i/Nx-0.5)^2 -(j/Ny-0.5)^2)/0.1);
    end
end
name = "measurement_data"+ string(Nx) + ".mat";
save(name, "psi")
end