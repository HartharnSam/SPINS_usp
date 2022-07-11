function qsp_to_physical(ii, var1, var2, spat_lims, var_lims, Region)
%QSP_TO_PHYSICAL - Relate region of qsp graph to physical space.
% Can be used interactively by clicking a region, or by setting it in the
% code
%
% Syntax:  qsp_to_physical(ii, var1, var2, spat_lims, var_lims, Region)
%
% Inputs:
%    ii - Simulation timestep to output for
%    var1 - variable to compare to var2 (on x axis)
%    var2 - Variable to compare to var1 (on y axis)
%    spat_lims - [optional] Spatial Region of physical space [xmin xmax zmin zmax]
%       optionally, only set the x limits. Defaults to full size of tank
%    var_lims - [optional] realistic limits of the variables to investigate as [var2min var2max var1min var2max]
%    Region - [optional] QSP Region of interest to display data for
%
% Other m-files required: qsp_mapped, spins_params, xgrid_reader,
% zgrid_reader, spins_reader_new, nearest_index, cmocean, plasma,
% figure_print_format
%
% See also: qsp_mapped
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 17-Jun-2022; Last revision: 17-Jun-2022
% MATLAB Version: 9.12.0.1956245 (R2022a) Update 2

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

% Option to invert the selected region (show outside the rectangle)
isInvert = false;
%% Load in data
% Compute the qsp data
params = spins_params;
if nargin<4
    xlims = [params.min_x params.min_x+params.Lx];
    spat_lims = xlims;
else
    xlims = spat_lims([1 2]);
end
if nargin >=4 && numel(spat_lims) == 4
    zlims = spat_lims([3 4]);
else
    zlims = [params.min_z params.min_z+params.Lz];
    spat_lims([3 4]) = zlims;
end

if nargin>4
    [qsp, myVar1, myKE, var_lims] = qsp_mapped(ii, var1, var2, spat_lims, var_lims);
else
    [qsp, myVar1, myKE, var_lims] = qsp_mapped(ii, var1, var2, spat_lims);
end

%% Set region of interest
if nargin > 5
    var1_ROI = Region([1 2]);
    var2_ROI = Region([3 4]);
else

    % Uncomment to set the region of interest interactively
    disp('Click the corners on the QSP plot to select your region of interest')
    [var1_ROI, var2_ROI] = ginput(2);
    var1_ROI = sort(var1_ROI); var2_ROI = sort(var2_ROI); % Sorting means the clicking can be in any order
end
%% Load in physical data
% read in the grids & cut down
x = xgrid_reader();
xmin_ind = nearest_index(x(:, 1), xlims(1));
xmax_ind = nearest_index(x(:, 1), xlims(2));
x = x(xmin_ind:xmax_ind, :);
z = zgrid_reader(xmin_ind:xmax_ind, []);

% note we still read in any z data from the region we want to cut due to
% the complicated cheb grid
switch var1
    case 's'
        data1 = spins_reader_new('s',ii, xmin_ind:xmax_ind, []);
        data1 = data1.*(data1>0); % We always know for salinity -ve isn't real
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
        u = spins_reader_new('u',ii, xmin_ind:xmax_ind, []);
        w = spins_reader_new('w',ii, xmin_ind:xmax_ind, []);
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
        u = spins_reader_new('u',ii, xmin_ind:xmax_ind, []);
        w = spins_reader_new('w',ii, xmin_ind:xmax_ind, []);
        data2 = sqrt(u.^2 + w.^2);
    otherwise
        try
            data2 = spins_reader_new(var2, ii, xmin_ind:xmax_ind, []);
        catch
            error([var2, ' not configured']);
        end
end

%% Extract QSP region of interest from physical data
RegOfInterest = (~((data1 >= var1_ROI(1)) & (data1 <= var1_ROI(2)) & (data2 >= var2_ROI(1))...
    & (data2 <= var2_ROI(2))));
if isInvert
    RegOfInterest = ~RegOfInterest;
end
data1(RegOfInterest) = NaN;
data2(RegOfInterest) = NaN;


%% Plot it up
% These are the upper plots of the variables in physical space
% Copy over from the figure made by QSP_mapped
hf1 = gcf;
hf2 = figure(2);
hf2.Position = hf1.Position;
ch = get(hf1, 'children');
nh = copyobj(ch, hf2);

ax1 = hf2.Children(end);
ax2 = hf2.Children(end-2);
ax3 = hf2.Children(end-3);
axes(ax1)
hold off
pcolor(ax1, x, z, data1); shading flat;
title(['t = ', num2str(ii)]);
colormap(gca, cmocean('dense'));
caxis(var1_ROI);
c = colorbar('location', 'EastOutside');
ylabel(c, var1); ylabel('z (m)');
axis([xlims zlims])
xticklabels([]);
hold on;
plot(x(:, 1), z(:, 1), 'k-');

axes(ax2);
hold off
pcolor(ax2, x, z, data2); shading flat;
if strcmpi(var2, 'vorty')
    colormap(gca, cmocean('balance'));
else
    colormap(gca, cmocean('amp'))
end
caxis(var2_ROI);
c = colorbar('Location', 'EastOutside');
ylabel(c, var2); xlabel('x (m)'); ylabel('z (m)');
hold on;
plot(x(:, 1), z(:, 1), 'k-');
axis([xlims zlims])
ax2.Position(3) = ax1.Position(3);

% Plot the QSP Region of Interest
axes(ax3);
hold on
rectangle('Position', [var1_ROI(1) var2_ROI(1) diff(var1_ROI) diff(var2_ROI)],...
    'EdgeColor', 'w');

%% Finish off the plot by formatting it all
figure_print_format;
