function ROI = usp_to_physical(ii, var1, var2, spatLims, varLims, region)
%USP_TO_PHYSICAL - Relate region of usp graph to physical space.
% Can be used interactively by clicking a region, or by setting it in the
% code
%
% Syntax:  usp_to_physical(ii, var1, var2, spatLims, varLims, region)
%
% Inputs:
%    ii - Simulation timestep to output for
%    var1 - variable to compare to var2 (on x axis)
%    var2 - Variable to compare to var1 (on y axis)
%    spatLims - [optional] Spatial region of physical space [xmin xmax zmin zmax]
%       optionally, only set the x limits. Defaults to full size of tank
%    varLims - [optional] realistic limits of the variables to investigate as [var2min var2max var1min var2max]
%    region - [optional] USP region of interest to display data for
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
% 17-Jun-2022; Last revision: 17-Jun-2022
% MATLAB Version: 9.12.0.1956245 (R2022a) Update 2
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

[~, ~, ~] = usp_2d(ii, var1, var2, spatLims, varLims, (nargout == 0));

%% Set region of interest
if nargin > 5
    var1ROI = region([1 2]);
    var2ROI = region([3 4]);
else
    
    % to set the region of interest interactively
    disp('Click the corners on the QSP plot to select your region of interest')
    [var1ROI, var2ROI] = ginput(2);
    var1ROI = sort(var1ROI); var2ROI = sort(var2ROI); % Sorting means the clicking can be in any order
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

%% Extract USP region of interest from physical data
RegOfInterest = (~((data1 >= var1ROI(1)) & (data1 <= var1ROI(2)) & (data2 >= var2ROI(1))...
    & (data2 <= var2ROI(2))));

if isInvert
    RegOfInterest = ~RegOfInterest;
end

data1(RegOfInterest) = NaN;
data2(RegOfInterest) = NaN;

if nargout == 0
    %% Plot it up
    % These are the upper plots of the variables in physical space
    % Copy over from the figure made by QSP_mapped
    hf1 = gcf;
    hf2 = figure(2);
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
    colormap(gca, cmocean('dense'));
    caxis(var1ROI);
    c = colorbar('location', 'EastOutside');
    ylabel(c, var1); ylabel('z (m)');
    axis([xlims zlims])
    xticklabels([]);
    hold on;
    plot(x(:, 1), z(:, 1), 'k-');
    
    % Replace the variable 2 plot with the ROI plots
    axes(ax2);
    hold off
    pcolor(ax2, x, z, data2); shading flat;
    if strcmpi(var2, 'vorty')
        colormap(gca, cmocean('balance'));
    else
        colormap(gca, cmocean('amp'))
    end
    caxis(var2ROI);
    c = colorbar('Location', 'EastOutside');
    ylabel(c, var2); xlabel('x (m)'); ylabel('z (m)');
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
    ROI.region = RegOfInterest;
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
