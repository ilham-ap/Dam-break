
%% Initial Program
clear;
clc;
clf;

%% Initial Condition
g       =9.81;                  % gravitasi
n       =0;                     % koefisien manning
Lx      =2;                   % panjang arah x
Ly      =0.16;                  % panjang arah y
px1     =0.9;                  % koordinat pilar arah-x
px2     =1.02;                  % koordinat pilar arah-x
py1     =0.3;                   % koordinat pilar arah-y
py2     =0.42;                  % koordinat pilar arah-y
fk      =0.99;
ct      =0.002;
dx      =0.005;                  % beda jarak arah x
dt      =ct*dx;                 % beda waktu
dy      =0.005;                  % beda jarak arah y
im      =Lx/dx+1;               % jumlah grid dalam arah x
jm      =Ly/dy+1;               % jumlah grid dalam arah y
time    =1;                   % lama simulasi
tm      =time/dt;               % jumlah loop waktu
tout    =0.001;                  % waktu simulasi
tplot   =tout/dt;               % waktu plot simulasi

% Mesh Grid
[x, y] = meshgrid(linspace(0,Ly,jm),linspace(0,Lx,im));

% Initial Value
h=zeros(im,jm);        un=zeros(im,jm);        vn=zeros(im,jm);
hn=zeros(im,jm);       u=zeros(im,jm);         v=zeros(im,jm);   
hfilter=zeros(im,jm);  ufilter=zeros(im,jm);   vfilter=zeros(im,jm); 
z=zeros(im,jm);

%% Initial Condition
for i=1:im
    for j=1:jm
        if i<0.4/dx+1
            h(i,j)=0.15;
        else 
            h(i,j)=0.015;
        end
    end
end

figure(1)
meshz(x,y,h);
colormap('gray')
shading('interp')
view(121,63)
xlim([0 Ly])
ylim([0 Lx])
zlim([0 0.4])
xlabel('X-axis (m)')
ylabel('Y-axis (m)')
zlabel('Z-axis (m)')
title('Initial Condition')


%% Main Program (FTCS)
for t = 1:tm+1
    ts = t*dt;
    for i = 2 : im
        for j = 2:jm
        % Kontinuitas
        hn(i,j)=h(i,j)-dt/dx*(u(i,j)*h(i,j)-u(i-1,j)*h(i-1,j))-dt/dy*(v(i,j)*h(i,j)-v(i,j-1)*h(i,j-1));
        % Momentum
        % arah x
        un(i,j)=u(i,j)-dt/dx*u(i,j)*(u(i,j)-u(i-1,j))-dt/dy*v(i,j)*(u(i,j)-u(i,j-1))-g*dt/dx*(h(i,j)-h(i-1,j));
        % arah y
        vn(i,j)=v(i,j)-dt/dy*v(i,j)*(v(i,j)-v(i,j-1))-dt/dx*u(i,j)*(v(i,j)-v(i-1,j))-g*dt/dy*(h(i,j)-h(i,j-1));
        end
    end
    
%%  Syarat Batas
% Syarat batas dinding
hn(:,1) = hn(:,2);      un(:,1) = un(:,2);      vn(:,1) = 0;
hn(:,jm) = hn(:,jm-1);  un(:,jm) = un(:,jm-1);  vn(:,jm) = 0;
hn(1,:) = hn(2,:);      un(1,:) = 0;            vn(1,:) = vn(2,:);
hn(im,:) = hn(im-1,:);  un(im,:) = 0;           vn(im,:) = vn(im-1,:);
    
      
%% Filter 2D
for i=2:im-1
    for j=2:jm-1
        hfilter(i,j)=fk*hn(i,j)+(1-fk)*(hn(i-1,j)+hn(i+1,j)+hn(i,j+1)+hn(i,j-1))/4;
        ufilter(i,j)=fk*un(i,j)+(1-fk)*(un(i-1,j)+un(i+1,j)+un(i,j+1)+un(i,j-1))/4;
        vfilter(i,j)=fk*vn(i,j)+(1-fk)*(vn(i-1,j)+vn(i+1,j)+vn(i,j+1)+vn(i,j-1))/4;
    end
end
    
% Syarat Batas Dinding
% Boundary dinding dinding
hfilter(:,1) = hfilter(:,2);      ufilter(:,1) = ufilter(:,2);      vfilter(:,1) = 0;
hfilter(:,jm) = hfilter(:,jm-1);  ufilter(:,jm) = ufilter(:,jm-1);  vfilter(:,jm) = 0;
hfilter(1,:) = hfilter(2,:);      ufilter(1,:) = 0;                 vfilter(1,:) = vfilter(2,:);
hfilter(im,:) = hfilter(im-1,:);  ufilter(im,:) = 0;                vfilter(im,:) = vfilter(im-1,:);

    
%% Loop Waktu
hn=hfilter;
un=ufilter;
vn=vfilter;
h=hn;
u=un;
v=vn;


if  mod(t,tplot)==0
    figure(2)
    meshz(x,y,h)
    title(['t = ' num2str(ts) ' second '])
    colormap('gray')
    shading('interp')
    view(121,63)
    xlim([0 Ly])
    ylim([0 Lx])
    zlim([0.015 0.4])
    xlabel('X-axis (m)')
    ylabel('Y-axis (m)')
    zlabel('Water Depth (m)')
    grid on
    pause(0.001)
    
    if ts==0.05
        figure(3)
    meshz(x,y,h)
    title(['t = ' num2str(ts) ' second '])
    colormap('gray')
    shading('interp')
    view(121,63)
    xlim([0 Ly])
    ylim([0 Lx])
    zlim([0 0.4])
    xlabel('X-axis (m)')
    ylabel('Y-axis (m)')
    zlabel('Water Depth (m)')
    grid on
    
    figure(4)
    contour(x,y,h)
    
    end
    
    if ts==0.1
        figure(5)
    meshz(x,y,h)
    title(['t = ' num2str(ts) ' second '])
    colormap('gray')
    shading('interp')
    view(121,63)
    xlim([0 Ly])
    ylim([0 Lx])
    zlim([0 0.4])
    xlabel('X-axis (m)')
    ylabel('Y-axis (m)')
    zlabel('Water Depth (m)')
    grid on
    
    figure(6)
    contour(x,y,h)
    
    end
    
    if ts==1
        figure(7)
    meshz(x,y,h)
    title(['t = ' num2str(ts) ' second '])
    colormap('gray')
    shading('interp')
    view(121,63)
    xlim([0 Ly])
    ylim([0 Lx])
    zlim([0 0.4])
    xlabel('X-axis (m)')
    ylabel('Y-axis (m)')
    zlabel('Water Depth (m)')
    grid on
    
    figure(8)
    contour(x,y,h)
    
    end
end
end