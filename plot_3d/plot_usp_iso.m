%% PLOT_USP_ISO - isosurface maker for USP tagged regions
% Plots 3D isosurfaces produced by make_usp_isos
% Loads in .mat file saved by the make_usp_isos script and plots them,
% for visualisation of 3D qsp_to_physical outputs
%

clearvars; close all; clc
%% Load data
% Manually create 3D unmapped grid 
Lx = 0.512; Ly = 0.128; Lz = 0.128;
Nx = 512;Ny = 128;Nz = 256;
x1d = (0.5:Nx-0.5)*Lx/Nx;
y1d = (0.5:Ny-0.5)*Ly/Ny;
z1d = (0.5:Nz-0.5)*Lz/Nz;
[x,y,z] = meshgrid(x1d,y1d,z1d);

% Load the .mat file & permute for Matlab's graphics ordering
load t_kev_ROI.mat
myreg = permute(RegOfInterest,[2,1,3]); % my region of interest

%% Plot data
figure(2)
clf
set(gcf,'DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
    'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
    'DefaultAxesFontWeight','bold');

%% Plot the isoline as reference
% Load in that data
rho = eqn_of_state(spins_reader_new('t', 45), 0);
rho = permute(rho, [2 1 3]);
% find max and min so you can better set the value to plot
mx=max(rho(:));
mn=min(rho(:));
delrho = mx-mn;

% Plot
p2b=patch(isosurface(x,y,z,rho,mn+(0.35*delrho)));
set(p2b,'FaceColor', [.5 .5 .5],'EdgeColor','none','FaceAlpha',0.7);
hold on

%% Plot the region of interest as isosurface
p2a = patch(isosurface(x,y,z,myreg,0.99));
set(p2a,'FaceColor','red','EdgeColor','none','FaceAlpha', .9);

%% Format the plot
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
%axis([0 2 0 6 0 20])
grid on
azimuth = -40; altitude = 12;
view(azimuth, altitude) % from high on left wall
%view(-60,32) % from the front
axis equal
lighting gouraud
camlight(-45,0.1) %top
xlim([x(1) x(end)]);
ylim([y(1) y(end)]);
zlim([z(1) z(end)]);

%% Make a funky video
% pause
% 
% for i = 0:2:360
%     view(azimuth+i, altitude);
%     pause(0.1);
% end