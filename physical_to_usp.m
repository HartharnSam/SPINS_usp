function physical_to_usp(ii, var1, var2, spatLims, varLims)
%USP_TO_PHYSICAL - Reverse idea of usp_to_physical - click a grid location
%and see where it sits in USP space, and identify other regions in the same
%USP space. Not sure of it's utility... yet...
% Used interactively by clicking a region
%
% Syntax:  physical_to_usp(ii, var1, var2, spatLims, varLims)
%
% Inputs:
%    ii - Simulation timestep to output for
%    var1 - variable to compare to var2 (on x axis)
%    var2 - Variable to compare to var1 (on y axis)
%    spatLims - [optional] Spatial region of physical space [xmin xmax zmin zmax]
%       optionally, only set the x limits. Defaults to full size of tank
%    varLims - [optional] realistic limits of the variables to investigate as [var2min var2max var1min var2max]
%
% Other m-files required: usp_2d, spins_params, xgrid_reader,
% zgrid_reader, spins_reader_new, nearest_index, cmocean, plasma,
% figure_print_format
%
% See also: usp_2d
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 23-Apr-2024; Last revision: 23-Apr-2024
% MATLAB Version: 9.10.0.1739362 (R2021a) Update 5
%
%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
%close all;
% Option to invert the selected region (show outside the rectangle)
isInvert = false;
%% Load in data
% Compute the qsp data
params = spins_params;
if nargin < 4 || isempty(spatLims)
    xlims = [params.min_x params.min_x+params.Lx];
    spatLims = xlims;
else
    xlims = spatLims([1 2]);
end
if nargin >= 4 && numel(spatLims) == 4
    zlims = spatLims([3 4]);
else
    zlims = [params.min_z params.min_z+params.Lz];
    spatLims([3 4]) = zlims;
end

if nargin <= 4
    varLims = [];
end

%% Load in physical data
% read in the grids & cut down
x = xgrid_reader();
xminInd = nearest_index(x(:, 1), xlims(1));
xmaxInd = nearest_index(x(:, 1), xlims(2));
x = x(xminInd:xmaxInd, :);
z = zgrid_reader(xminInd:xmaxInd, []);

% note we still read in any z data from the region we want to cut due to
% the complicated cheb grid
switch var1
    case 's'
        data1 = spins_reader_new('s',ii, xminInd:xmaxInd, []);
        data1 = data1.*(data1 > 0); % We always know for salinity -ve isn't real
    case 'rho'
        data1 = spins_reader_new('rho', ii, xminInd:xmaxInd, []);
    case 'rho_z2'
        try
            data1 = spins_reader_new('rho_z', ii, xminInd:xmaxInd, []).^2;
        catch
            spins_derivs('rho_z', ii, true);
            data1 = spins_reader_new('rho_z', ii, xminInd:xmaxInd, []).^2;
        end
    otherwise
        try
            data1 = spins_reader_new(var1, ii, xminInd:xmaxInd, []);
        catch
            error([var1, ' not configured']);
        end
end

switch lower(var2)
    case 'ke'
        u = spins_reader_new('u',ii, xminInd:xmaxInd, []);
        w = spins_reader_new('w',ii, xminInd:xmaxInd, []);
        data2 = 0.5*(u.^2+w.^2);
    case 'vorty'
        data2 = spins_reader_new('vorty', ii, xminInd:xmaxInd, []);
    case 'enstrophy'
        try
            data2 = spins_reader_new('enstrophy', ii, xminInd:xmaxInd, []);
        catch
            data2 = 0.5*spins_reader_new('vorty', ii, xminInd:xmaxInd, []).^2;
        end
    case 'diss'
        data2 = spins_reader_new('diss', ii, xminInd:xmaxInd, []);
        data2 = log10(data2);
    case 'speed'
        u = spins_reader_new('u',ii, xminInd:xmaxInd, []);
        w = spins_reader_new('w',ii, xminInd:xmaxInd, []);
        data2 = sqrt(u.^2 + w.^2);
    otherwise
        try
            data2 = spins_reader_new(var2, ii, xminInd:xmaxInd, []);
        catch
            error([var2, ' not configured']);
        end
end

%% And make the original plot
tiledlayout(2, 1);
ax1 = nexttile; ax2 = nexttile;

pcolor(ax1, x, z, data1); colorbar(ax1);
shading(ax1, 'flat');

title(ax1, ['t = ', num2str(ii)]);
colormap(ax1, cmocean('dense'));
c = colorbar(ax1, 'location', 'EastOutside');
ylabel(c, var1); ylabel(ax1, 'z (m)');
axis(ax1, [xlims zlims])
xticklabels(ax1, []);
hold(ax1, 'on');
plot(ax1, x(:, 1), z(:, 1), 'k-');

pcolor(ax2, x, z, data2); colorbar(ax2);
hold(ax2,'off');
shading(ax2, 'flat');
if strcmpi(var2, 'vorty')
    colormap(ax2, cmocean('balance'));
else
    colormap(ax2, cmocean('amp'))
end
c = colorbar(ax2, 'Location', 'EastOutside');
ylabel(c, var2); xlabel(ax2, 'x (m)'); ylabel(ax2, 'z (m)');
hold(ax2, 'on');
plot(ax2, x(:, 1), z(:, 1), 'k-');
axis(ax2, [xlims zlims])


%% Set region of interest
axes(ax1);
set(gcf, 'Position', groot().MonitorPositions(end, :));
disp('Click the corners on the QSP plot to select your region of interest')
[xROI, zROI] = ginput(1);

xROI_ind = nearest_index(x(:, 1), xROI);
zROI_ind = nearest_index(z(xROI_ind, :), zROI);

data1ROI = data1(xROI_ind, zROI_ind);
data2ROI = data2(xROI_ind, zROI_ind);

%% Plot the picked data on a standard USP plot

clf;
[~, myVar1, myVar2] = usp_2d(ii, var1, var2, spatLims, varLims, true);

hf1 = gcf;
aces = findobj(hf1,'Type','Axes');
[~, aces] = sort_axes(aces);

ax1USP = aces(1);
ax2USP = aces(2);
ax3_usp = aces(4);
hold(ax3_usp, 'on')
plot(ax3_usp, data1ROI, data2ROI, 'xw');
plot(ax1USP, xROI, zROI, 'xw');
plot(ax2USP, xROI, zROI, 'xw');

%% plot a ROI plot based on the single bin
% but first, identify the bin
data1Bin = interp1(myVar1, 1:length(myVar1), data1ROI);
data1Bin = myVar1(floor(data1Bin) + 0:1);

data2Bin = interp1(myVar2, 1:length(myVar2), data2ROI);
data2Bin = myVar1(floor(data2Bin) + 0:1);

% Run usp_to_physical, feeding in the single-bin ROI
figure
usp_to_physical(ii, var1, var2, spatLims, varLims, [data1Bin data2Bin]);

end

function [sortedPositions, sortedAxes] = sort_axes(arrayOfAxes)
% SORT_AXES sorts the axis from top-left to bottom-right.
% [POSITIONS,AXES] = sort_axes(arrayOfAxes) Takes in an array of subplot axes
% and sorts them from top-left to bottom right according to their position.
% This returns POSITIONS which is a matrix that contains the position
% vectors of the sorted axes. AXES is the array of sorted axes.
numAxes = length(arrayOfAxes);
positions = zeros(numAxes,4);
for ii = 1:numAxes
    positions(ii,:) = arrayOfAxes(ii).Position;
end
[sortedPositions,sortIndex] = sortrows(positions,[-2 1]);
sortedAxes = arrayOfAxes(sortIndex);
end
