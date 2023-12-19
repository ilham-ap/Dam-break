clc;
clear;

Lx=2;
Ly=2;
nx=41;
ny=41;
dx=Lx/(nx-1);
dy=Ly/(ny-1);
g=9.81;
[x y] = meshgrid(linspace(0,Lx,nx),linspace(0,Ly,ny));

eta=zeros(nx,ny);  up=zeros(nx,ny);  vp=zeros(nx,ny);
etam=zeros(nx,ny);  u=zeros(nx,ny);   v=zeros(nx,ny);   
etap=zeros(nx,ny); um=zeros(nx,ny);  vm=zeros(nx,ny);

eta(:,1) = eta(:,2);      up(:,1) = up(:,2);       vp(:,1) = vp(:,2);
eta(:,ny) = eta(:,ny-1);  up(:,ny) = up(:,ny-1);   vp(:,ny) = vp(:,ny-1);
eta(1,:) = eta(2,:);      up(1,:) = up(2,:);      vp(1,:) = vp(2,:);
eta(nx,:) = eta(nx-1,:);  up(nx,:) = up(nx-1,:);  vp(nx,:) = vp(nx-1,:);

for i=1:nx 
   for j=1:ny
   if y(i,j) > 0.4;
       etam(i,j)=0;
   else
       etam(i,j)=0.15;
   end
       
   end
end
eta=etam;
dt=.001;
Nsteps=100;
for n=1:Nsteps
    t=n*dt;
    
    for i=2:nx-1
        for j=2:ny-1
    up(i,j)=um(i,j)-(dt/dx)*um(i,j)*(um(i,j)-um(i-1,j))-(dt/dy)*vm(i,j)*(um(i,j)-um(i,j-1))-g*(dt/dx)*(eta(i,j)-eta(i-1,j));
    vp(i,j)=vm(i,j)-(dt/dy)*vm(i,j)*(vm(i,j)-vm(i,j-1))-(dt/dx)*um(i,j)*(vm(i,j)-vm(i-1,j))-g*(dt/dy)*(eta(i,j)-eta(i,j-1));
    %etap(i,j)=eta(i,j)-(dt/dx)*(eta(i,j)-eta(i-1,j))*(um(i,j)-um(i-1,j))-(dt/dy)*(eta(i,j)-eta(i,j-1)*(vm(i,j)-vm(i,j-1))); 
    
   etap(i,j)=eta(i,j)-(dt/dx)*(eta(i,j)-eta(i-1,j))-(dt/dy)*(eta(i,j)-eta(i,j-1));
    
        end
    end
    um=up;
    vm=vp;
 eta=etap;
figure(1)
z=etap;
surf(x,y,etap)
colormap('Gray')
zlim([-0.4 0.4])
shading interp
zlabel('H +/- \eta')
xlabel('x')
ylabel('y')
end