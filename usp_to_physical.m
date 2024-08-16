function ROI = usp_to_physical(ii, var1, var2, spatLims, varLims, region, opts)
%USP_TO_PHYSICAL - Maps a selected region from a USP graph to physical space.
%   Can be used interactively by clicking a region in the USP graph, or
%   programmatically by specifying the region in the arguments
%
% Syntax:
%   ROI = usp_to_physical(ii, var1, var2, spatLims, varLims, region, opts)
%
% Inputs:
%    ii         - Index of the simulation output to process
%    var1       - First variable to compare (plotted on x axis)
%    var2       - Second variable to compare (plotted on y axis)
%    spatLims   - [optional] Spatial limits in the format [xmin xmax zmin zmax]
%               Defaults to full size of tank, optionally only [xmin xmax]
%    varLims    - [optional] limits of variables to investigate in the format
%               [var2min var2max var1min var2max]. Note: the order is var2
%               (y-axis), var1 (x-axis) as it is more common to set only
%               the y limits [var2min var2max] than only x!
%    region     - [optional] USP region of interest to display data for. If
%               empty, the user will be prompted to interactively provide
%               it
%    opts       - [Optional] Variable-Value pairs to input data directly
%
%    ROI - Region of Interest in Physical space corresponding to selected USP region
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
% 17-Jun-2022; Last revision: 17-Aug-2024
% MATLAB Version: 9.12.0.1956245 (R2022a) Update 2
%
%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
arguments
    ii (1, 1) uint16
    var1 (1, 1) string
    var2 (1, 1) string
    spatLims (:, 1) double = []
    varLims (:, 1) double = []
    region (:, 1) double = []
    opts.data1 (:, :) double = []
    opts.data2 (:, :) double = []
    opts.isInvert (1, 1) logical = false % Option to invert the selected region (show outside the rectangle)
end

%% Load in grid data
params = spins_params;
if isempty(spatLims)
    xlims = [0 params.Lx]+params.min_x;
    spatLims = xlims;
else
    xlims = spatLims([1 2]);
end
if numel(spatLims) == 4
    zlims = spatLims([3 4]);
else
    zlims = [0 params.Lz]+params.min_z;
    spatLims([3 4]) = zlims;
end

% read in the grids & cut down
x = xgrid_reader();
xminInd = nearest_index(x(:, 1), xlims(1));
xmaxInd = nearest_index(x(:, 1), xlims(2));
z = zgrid_reader(xminInd:xmaxInd, []);

if isequal(params.mapped_grid, 'true')
    zInds = [];
    x = x(xminInd:xmaxInd, :);
    z = z(xminInd:xmaxInd, :);
else
    zminInd = nearest_index(z(1, :), zlims(1));
    zmaxInd = nearest_index(z(1, :), zlims(2));
    zInds = zminInd:zmaxInd;
    x = x(xminInd:xmaxInd, zInds);
    z = z(xminInd:xmaxInd, zInds);
end

%% Load in physical data

if isempty(opts.data1)
    data1 = get_spins_data(var1, ii, xminInd, xmaxInd, zInds);
else
    if isequal(size(x), size(opts.data1))
        data1 = opts.data1;
    else
        data1 = opts.data1(xminInd:xmaxInd, zInds);
    end
end
if isempty(opts.data2)
    data2 = get_spins_data(var2, ii, xminInd, xmaxInd, zInds);
else
    if isequal(size(x), size(opts.data2))
        data2 = opts.data2;
    else
        data2 = opts.data2(xminInd:xmaxInd, zInds);
    end
end

[~, ~, ~] = usp_2d(ii, var1, var2, spatLims, varLims, (nargout == 0), "data1", data1, "data2", data2);

%% Set region of interest
if (nargin > 5)
    var1ROI = region([1 2]);
    var2ROI = region([3 4]);
else
    set(gcf, 'Position', groot().MonitorPositions(end, :));
    disp('Click the corners on the QSP plot to select your region of interest')
    [var1ROI, var2ROI] = ginput(2);
    var1ROI = sort(var1ROI); var2ROI = sort(var2ROI); % Sorting means the clicking can be in any order
end

%% Extract USP region of interest from physical data
regionOfInterest = (~((data1 >= var1ROI(1)) & (data1 <= var1ROI(2)) & (data2 >= var2ROI(1))...
    & (data2 <= var2ROI(2))));

if opts.isInvert
    regionOfInterest = ~regionOfInterest;
end

data1(regionOfInterest) = NaN;
data2(regionOfInterest) = NaN;

[axLab1] = get_axis_labels(var1);
[axLab2] = get_axis_labels(var2);

if (nargout == 0)
    %% Plot it up
    % These are the upper plots of the variables in physical space
    % Copy over from the figure made by QSP_mapped
    hf1 = gcf;
    hf2 = figure(hf1.Number +1);
    hf2.Position = hf1.Position;
    ch = get(hf1, 'children');
    copyobj(ch, hf2);
    %%
    aces = findobj(hf2,'Type','Axes');
    [~, aces] = sort_axes(aces);
    %%
    ax1 = aces(1);
    ax2 = aces(2);
    ax3 = aces(4);
    
    % Replace the variable 1 plot with the ROI plots
    axes(ax1)
    hold off
    pcolor(ax1, x, z, data1); shading flat;
    title(['t = ', num2str(ii)]);
    caxis(var1ROI);
    colormap(gca, cmocean('dense'));
    c = colorbar('location', 'EastOutside');
    ylabel(c, axLab1); ylabel('z (m)');
    axis([xlims zlims])
    xticklabels([]);
    hold on;
    plot(x(:, 1), z(:, 1), 'k-');
    
    % Replace the variable 2 plot with the ROI plots
    axes(ax2);
    hold off
    pcolor(ax2, x, z, data2); shading flat;
    caxis(var2ROI);
    if strcmpi(var2, 'vorty')
        colormap(gca, cmocean('balance'));
    else
        colormap(gca, cmocean('amp'))
    end
    c = colorbar('Location', 'EastOutside');
    ylabel(c, axLab2); xlabel('x (m)'); ylabel('z (m)');
    hold on;
    plot(x(:, 1), z(:, 1), 'k-');
    axis([xlims zlims])
    ax2.Position(3) = ax1.Position(3);
    
    % Plot the QSP region of Interest
    axes(ax3);
    hold on
    rectangle(ax3, 'Position', [var1ROI(1) var2ROI(1) diff(var1ROI) diff(var2ROI)],...
        'EdgeColor', 'w');
    
    %% Finish off the plot by formatting it all
    figure_print_format;
    
    %%
else
    ROI.region = regionOfInterest;
    ROI.x = x;
    ROI.z = z;
end
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
