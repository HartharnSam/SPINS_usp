%% PLOT_QSP_ISO - Marek's isosurface maker for Sam's QSP tagged regions
% Plots 3D isosurfaces produced by make_qsp_isos
% Loads in .mat file saved by the make_qsp_isos script and plots them,
% for visualisation of 3D qsp_to_physical outputs



%% Read and permute for Matlab's graphics ordering
%load ROI_green_45s.mat
load ROI_white_45s.mat

Lx=0.512,Ly=0.128;Lz=0.128;
Nx=512;Ny=128;Nz=256;
x1d=(0.5:Nx-0.5)*Lx/Nx;
y1d=(0.5:Ny-0.5)*Ly/Ny;
z1d=(0.5:Nz-0.5)*Lz/Nz;

%t=spins_reader_new('t'45);
%rho=eqn_of_state(t,0*t);

%load mydenfor3d_norot.mat

%rho=permute(denfor3d,[2,1,3]);
myreg=permute(RegOfInterest,[2,1,3]);

clear t

% find max and min so you can better set the value to plot
%mx=max(rho(:));
%mn=min(rho(:));
% grid
[x,y,z]=meshgrid(x1d,y1d,z1d);



%An example with multiple isosurface values
figure(2)
clf
set(gcf,'DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
    'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
     'DefaultAxesFontWeight','bold');
%subplot(1,2,1) %for report Figure
%p2a=patch(isosurface(x,y,z,rho,0.4*(mn+mx)));
p2a=patch(isosurface(x,y,z,myreg,0.99));

set(p2a,'FaceColor','red','EdgeColor','none','FaceAlpha',0.3);
 hold on
% p2b=patch(isosurface(x,y,z,rho,0.45*(mn+mx)));
%  p2b=patch(isosurface(x,y,z,rho,999.972));
%  set(p2b,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.7);
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
%axis([0 2 0 6 0 20])
grid on
view(40,52) % from high on left wall
view(-60,32) % from the front

lighting gouraud
camlight(-45,0.1) %top
camlight(0,1) %top