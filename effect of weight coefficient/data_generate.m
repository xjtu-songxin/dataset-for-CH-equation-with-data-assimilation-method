close all;clc;clear 
rng('default') 
rng(1);
Nx = 100;
Ny = Nx;
Nz = Nx;
xa = 0; xb = 1;
ya = xa ; yb = xb;
za = xa ; zb = xb;
h=(xb-xa)/Nx;
x = linspace(xa+0.5*h,xb-0.5*h,Nx);
y = linspace(ya+0.5*h,yb-0.5*h,Ny);
z = linspace(za+0.5*h,zb-0.5*h,Nz);
dt = 1e-4;    
epsilon = 0.01;
maxiter = 5000;
ns= 500;

alpha_list= [1,0.01,0.05,0.1,0.5];

lambda = 50;
beta = 1;


[yy,xx,zz] = meshgrid(x,y,z);
kx = 2*pi/(xb-xa)*[0:Nx/2 -Nx/2+1:-1];
ky = 2*pi/(yb-ya)*[0:Ny/2 -Ny/2+1:-1];
kz = 2*pi/(zb-za)*[0:Nz/2 -Nz/2+1:-1];
k2x = kx.^2;   k2y = ky.^2; k2z = kz.^2;
k4x = kx.^4;   k4y = ky.^4; k4z = kz.^4;
[kyy,kxx,kzz]=meshgrid(k2x,k2y,k2z);
[kyyyy,kxxxx,kzzzz]=meshgrid(k4x,k4y,k4z);
wei=(kxx+kyy+kzz);
wei1=(kxxxx + kyyyy + kzzzz + 2.*kxx.*kyy + 2.*kxx.*kzz + 2.*kzz.*kyy );
 
for alpha =alpha_list


    
ave=0.4;
phi=1*(rand(Nx,Ny,Nz));
phi(phi>ave)= 1;
phi(phi<ave)=-1;

start_mass=sum(phi(:))/(Nx*Ny*Nz);


rng(2);
mask = rand(Nx,Ny,Nz);
theshold = 0;
mask(mask>theshold) = 1;
mask(mask<=theshold) = 0;

%% load measurement value
psi = zeros(Nx,Ny,Nz,1);
for ix=1:1:Nx
    for iy = 1:1:Ny
        for iz = 1:1:Nz
            psi(ix,iy,iz,1) = tanh(( (ix/Nx-0.5)^2 +(iy/Ny-0.5)^2+(iz/Nz-0.5)^2 -0.2025 )/0.1);
        end
    end
end
                
            
name = "measurement_data"+ string(Nx) + ".mat";
save(name, "psi")
end