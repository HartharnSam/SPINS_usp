function [qsp, myVar1, myVar2, var2_lims] = qsp_mapped(ii, var1, var2, xlims, var2_lims)
%FUNCTION_NAME - One line description of what the function or script performs (H1 line)%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [qsp, myVar1, myke, ke_lims] = qsp_mapped(ii, var1, xlims, var2_lims)
%
% Inputs:
%    ii - Simulation timestep to output for
%    var1 - variable to compare to var2 (on x axis)
%    var2 - Variable to compare to var1 (on y axis)
%    xlims - Region of tank to compare for
%    var2_lims - limits of variable 2 to investigate as [min max]
%
% Outputs:
%    qsp - Histogram of combined volume
%    myVar1 - variable values for each histogram box
%    myVar2 - ""
%    var2_lims - Limits of the second variable calculated by the script
%
% Other m-files required: xgrid_reader, zgrid_reader, spins_params,
% spins_reader_new, nearest_index, cheb, figure_print_format,
% plasma, subaxis
% Subfunctions: none
% MAT-files required: none
%
% See also: qsp_to_physical, spins_QSP_csv
% Author:
% Department of Applied Mathematics, University of Waterloo
% email address:
% GitHub:
% 14-Jun-2022; Last revision: 20-Jun-2022
% MATLAB Version: 9.12.0.1956245 (R2022a) Update 2

close all;

% read in the grids & cut down to size
params = spins_params;
if nargin<4
    xlims = [params.min_x params.min_x+params.Lx];
end

x = xgrid_reader();
xmin_ind = nearest_index(x(:, 1), xlims(1));
xmax_ind = nearest_index(x(:, 1), xlims(2));
x = x(xmin_ind:xmax_ind, :);
z = zgrid_reader(xmin_ind:xmax_ind, []);

% get the size of the grids
[Nx, Nz] = size(x);

% parameter for the chebyshev grid which are on [-1,1]
Nzc = Nz-1;
% read in data
switch lower(var1)
    case 's'
        data1 = spins_reader_new('s',ii, xmin_ind:xmax_ind, []);
        data1 = data1.*(data1>0);
    case 'rho'
        data1 = spins_reader_new('rho', ii, xmin_ind:xmax_ind, []);
    otherwise
        try
            data1 = spins_reader_new(var1, ii, xmin_ind:xmax_ind, []);
        catch
            error([var1, ' not configured']);
        end
end

switch lower(var2)
    case 'ke'
        u = spins_reader_new('u',ii, xmin_ind:xmax_ind, []); w = spins_reader_new('w',ii, xmin_ind:xmax_ind, []);
        data2 = 0.5*(u.^2+w.^2);
    case 'vorty'
        data2 = spins_reader_new('vorty', ii, xmin_ind:xmax_ind, []);
    case 'enstrophy'
        try
            data2 = spins_reader_new('enstrophy', ii, xmin_ind:xmax_ind, []);
        catch
            data2 = 0.5*spins_reader_new('vorty', ii, xmin_ind:xmax_ind, []).^2;
        end
    case 'diss'
        data2 = spins_reader_new('diss', ii, xmin_ind:xmax_ind, []);
        data2 = log10(data2);
    case 'speed'
        u = spins_reader_new('u',ii, xmin_ind:xmax_ind, []); w = spins_reader_new('w',ii, xmin_ind:xmax_ind, []);
        data2 = sqrt(u.^2 + w.^2);
    otherwise
        try
            data1 = spins_reader_new(var2, ii, xmin_ind:xmax_ind, []);
        catch
            error([var2, ' not configured']);
        end
end

%% Compute the area
% Compute the area associated with each Chebyshev point using the values
% halfway between the point below and above
[~, z1dc] = cheb(Nzc);

% first do it on the standard interval
% the bottom and top most pts get a half grid box
arc(1) = 0.5*(z1dc(1)-z1dc(2));
arc(Nzc+1) = arc(1);

% over the interior pts and store the grid boxes
arc(2:Nzc) = 0.5*(z1dc(1:end-2)-z1dc(2:end-1))+0.5*(z1dc(2:end-1)-z1dc(3:end));

% now for each x point, stretch or shrink according to the local depth
Lznow = max(z, [], 2) - min(z, [], 2);
arcphys = arc.*Lznow;

% and get a long vector of the areas and the total area
arcphysv = arcphys(:);
totar = sum(arcphysv);

%% Sort data into the histogram "boxes"
numpts = 50;
% for variable on x
if strcmpi(var1, 'rho')
    smin = -params.delta_rho/2;
    smax = -smin;
else
    smin = min(data1(:));
    smax = max(data1(:));
end
data1(data1<smin) = smin;
data1(data1>smax) = smax;

ranges = smax-smin;
ds = ranges/(numpts-1);
myVar1 = smin+(0.5:numpts-0.5)*ds;

% For variable on y (KE)
if nargin<=4
    kemin = min(data2(:));
    kemax = max(data2(:));
else
    kemin = var2_lims(1);
    kemax = var2_lims(2);
end
data2(data2<kemin) = kemin;
data2(data2>kemax) = kemax;

rangeke = kemax-kemin;
dke = rangeke/(numpts-1);
myVar2 = kemin+(0.5:numpts-0.5)'*dke;

myhist = zeros(numpts,numpts);

%figure out which box coordinate you are in
sbox = ceil((data1-smin)/ds);
sbox = sbox+1*(data1==smin);

kebox = ceil((data2-kemin)/dke);
kebox = kebox+1*(data2==kemin);

% brutally inefficient but will work for 2D
% double loop
for i = 1:Nx
    for jj = 1:Nz
        % update the corect box's total with the current area value
        myhist(sbox(i, jj),kebox(i, jj)) = myhist(sbox(i, jj), kebox(i, jj))+1*arcphys(i,jj);
    end
end
qsp = myhist'/totar;

%% Plot up
isSanityCheck = true;
if isSanityCheck % These are the upper plots of the variables in physical space
    figure(1)
    ax1 = subaxis(4, 1, 1, 'MT', .04);    
    pcolor(x, z, data1), shading flat;
    title(['t = ', num2str(ii)]);
    colormap(gca, cmocean('dense'));
    caxis([smin smax]);
    c = colorbar('Location','EastOutside');
    ylabel(c, var1); ylabel('z (m)');
    axis tight
    xticklabels([]);
    hold on;
    plot(x(:, 1), z(:, 1), 'k-');

    ax2 = subaxis(4, 1, 2, 'MT', 0.04);
    pcolor(x,z,data2), shading flat;
    if strcmpi(var2, 'vorty')
        colormap(gca, cmocean('balance'));
    else
        colormap(gca, cmocean('amp'))
    end
    caxis([kemin kemax]);
    c = colorbar('Location','EastOutside');
    ylabel(c, var2); xlabel('x (m)'); ylabel('z (m)');
    hold on;
    plot(x(:, 1), z(:, 1), 'k-');
    axis tight
else
    tiledlayout(2, 1);
end

% Change the figure aspect ratio to taller
fig = gcf; fig.Position([3 4]) = [643.2000 531.2000];

ax3 = subaxis(6, 3, 1, 4, 1, 2, 'MB', 0.04);
ax4 = subaxis(6, 3, 2, 4, 1, 2, 'MB', 0.04);
ax5 = subaxis(6, 3, 2, 6, 1, 1, 'MB', 0.04);

% Plot the 2d QSP histogram
axes(ax4)
pcolor(ax4, myVar1, myVar2, log10(qsp));
shading flat
yticklabels([]);
xticklabels([]);
colormap(gca, plasma);
caxis([-6 -2]);
c = colorbar; ylabel(c, 'Volume');
axis tight;
box on
ax4.Position = ax4.Position;

% 1D Histogram for variable 2
axes(ax3)
plot(ax3, sum(qsp, 2)./sum(qsp(:)), myVar2, 'r-');
axis tight;
xlim([0 .15]);
ylabel(var2);
ax3.Position = [ax4.Position(1)-ax5.Position(4) ax4.Position(2) ax5.Position(4) ax4.Position(4)];
xticklabels([]);

% 1D Histogram for variable 1
axes(ax5)
plot(myVar1, sum(qsp)./sum(qsp(:)), 'b-');
axis tight;
ylim([0 .15]);
xlabel(var1);
ax5.Position = [ax4.Position(1) ax4.Position(2)-ax5.Position(4) ax4.Position(3) ax5.Position(4)];
yticklabels([]);

if nargout > 3
    var2_lims = [kemin kemax];
end

figure_print_format(gcf);
