function physical_to_usp(ii, var1, var2, spatLims, varLims, opts)
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
arguments 
    ii (1, 1) uint16
    var1 (1, 1) string
    var2 (1, 1) string
    spatLims (1, :) double = []
    varLims (1, :) double = []
    opts.data1 (:, :) double = []
    opts.data2 (:, :) double = []
    opts.isInvert (1, 1) logical = false;% Option to invert the selected region (show outside the rectangle)
end

%% Load in data
% Compute the qsp data
params = spins_params;
if isempty(spatLims)
    xlims = [params.min_x params.min_x+params.Lx];
    spatLims = xlims;
else
    xlims = spatLims([1 2]);
end
if (numel(spatLims) == 4)
    zlims = spatLims([3 4]);
else
    zlims = [params.min_z params.min_z+params.Lz];
    %spatLims([3 4]) = zlims;
end

if nargin <= 4
    varLims = [];
end

%% Load in physical data
% read in the grids & cut down
x = xgrid_reader();
xminInd = nearest_index(x(:, 1), xlims(1));
xmaxInd = nearest_index(x(:, 1), xlims(2));
z = zgrid_reader(xminInd:xmaxInd, []);

isCheb = isequal(params.mapped_grid, 'true') || isequal(params.type_z, 'NO_SLIP');

if isequal(params.mapped_grid, 'true')
    zInds = [];
    x = x(xminInd:xmaxInd, :);
    %z = z(xminInd:xmaxInd, :);
    
else
    zminInd = nearest_index(z(1, :), zlims(1));
    zmaxInd = nearest_index(z(1, :), zlims(2));
    zInds = zminInd:zmaxInd;
    x = x(xminInd:xmaxInd, zInds);
    z = z(:, zInds);
end

%% Read in data
if isempty(opts.data1)
    data1 = get_spins_data(var1, ii, xminInd, xmaxInd, zInds);
else
    data1 = opts.data1;
end
if isempty(opts.data2)
    data2 = get_spins_data(var2, ii, xminInd, xmaxInd, zInds);
else
    data2 = opts.data2;
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
disp('Click the point of interest on the upper plot')
[xROI, zROI] = ginput(1);

xROI_ind = nearest_index(x(:, 1), xROI);
zROI_ind = nearest_index(z(xROI_ind, :), zROI);

data1(data1>varLims(4)) = varLims(4);
data2(data2>varLims(2)) = varLims(2);

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
plot(ax3_usp, data1ROI, data2ROI, 'xw', 'MarkerSize',10);
plot(ax3_usp, data1ROI, data2ROI, 'xk', 'MarkerSize', 5);

plot(ax1USP, xROI, zROI, 'xw');
plot(ax2USP, xROI, zROI, 'xw');

%% plot a ROI plot based on the single bin
% but first, identify the bin
data1Bin = interp1(myVar1, 1:length(myVar1), data1ROI);
data1Bin = myVar1(floor(data1Bin) + (0:1));

data2Bin = interp1(myVar2, 1:length(myVar2), data2ROI);
data2Bin = myVar2(floor(data2Bin) + (0:1))';

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
