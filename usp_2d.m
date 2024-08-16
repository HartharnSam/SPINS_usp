function [usp, myVar1, myVar2, varLims] = usp_2d(ii, var1, var2, spatLims, varLims, isPlot, opts)
%USP_2D - Generates a USP plot, or paired histogram to identify the spatial
% overlap and interaction between two variables for 2D data.
%
% Syntax:
%    [usp, myVar1, myVar2, varLims] = usp_2d(ti, var1, var2, spatLims,
%    varLims, isPlot, "data1", data1, "data2", data2);
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
%   isPlot      - [Optional] Boolean flag for making a plot (true by default)
%   opts        - [Optional] Variable-Value pairs to input data directly
%
% Outputs:
%    usp        - 2D histogram of combined volume
%    myVar1     - Values of var1 in each histogram bin
%    myVar2     - ""
%    varLims    - Limits of the second variable calculated in the function
%
% Example:
%   usp = usp_2d(15, 'rho', 'KE', [0 1 0 0.3], [], true, "data1", rho);
%
% Other m-files required: xgrid_reader, zgrid_reader, spins_params,
% spins_reader_new, nearest_index, cheb, figure_print_format,
% plasma, subaxis
% Subfunctions: none
% MAT-files required: none
%
% See also: usp_to_physical, spins_QSP_csv
% Author: Sam Hartharn-Evans
% Department of Geography & Environmental Sciences, Northumbria University
% email address: sam.hartharn-evans@northumbria.ac.uk
% GitHub: https://github.com/HartharnSam
% 14-Jun-2022; Last revision: 16-Aug-2024
% MATLAB Version: 9.12.0.1956245 (R2022a) Update 2
%
% Dev notes:
% - I think this won't work for anything that isn't FOURIER in x, cheb in
% z. Relatively minor changes needed if they are simply the other way
% around, more if they are both fourier?
%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
arguments
    ii (1, 1) uint16
    var1 (1, 1) string
    var2 (1, 1) string
    spatLims (1, :) double = []
    varLims (1, :) double = []
    isPlot (1, 1) logical = true;
    opts.data1 (:, :) double = []
    opts.data2 (:, :) double = []
end
%% read in the grids
% & cut down to size
params = spins_params;
if isempty(spatLims)
    xlims = [0 params.Lx]+params.min_x;
else
    xlims = spatLims([1 2]);
end
if (numel(spatLims) == 4)
    zlims = spatLims([3 4]);
else
    zlims = [0 params.Lz]+params.min_z;
end

% read in the grids & cut down
x = xgrid_reader();
xminInd = nearest_index(x(:, 1), xlims(1));
xmaxInd = nearest_index(x(:, 1), xlims(2));
z = zgrid_reader(xminInd:xmaxInd, []);

isCheb = isequal(params.mapped_grid, 'true') || isequal(params.type_z, 'NO_SLIP');

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
% get the size of the grids
[Nx, Nz] = size(x);

% parameter for the chebyshev grid which are on [-1,1]
Nzc = Nz-1;
%% read in data
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
%% Sort data into the histogram "boxes"
numpts = 50;
% for variable on x
if (nargin > 4) && (numel(varLims) == 4)
    var1min = varLims(3); var1max = varLims(4);
else % Use defaults for the var1 limits
    %if strcmpi(var1, 'rho')
    %    var1min = -params.delta_rho/2;
    %    var1max = -var1min;
    %else
    var1min = min(data1(:));
    var1max = max(data1(:));
    %end
end
% Cap the values at the limits - ignoring these values by NaN'ing them
% would be a reasonable choice, but requires some more coding to make it
% work
data1(data1 < var1min) = var1min;
data1(data1 > var1max) = var1max;

ranges = var1max-var1min;
dVar1 = ranges/(numpts-1);
myVar1 = var1min+(0.5:numpts-0.5)*dVar1;

% For variable on y (KE)
if (nargin <= 4) || (isempty(varLims)) % Use defaults
    var2min = min(data2(:));
    var2max = max(data2(:));
else
    var2min = varLims(1);
    var2max = varLims(2);
end
data2(data2 < var2min) = var2min;
data2(data2 > var2max) = var2max;

rangeVar2 = var2max-var2min;
dVar2 = rangeVar2/(numpts-1);
myVar2 = var2min+(0.5:numpts-0.5)'*dVar2;

myhist = zeros(numpts,numpts);

%figure out which box coordinate you are in
var1box = ceil((data1-var1min)/dVar1);
var1box = var1box+1*(data1 == var1min);

var2box = ceil((data2-var2min)/dVar2);
var2box = var2box+1*(data2 == var2min);

if isCheb
    %% Compute the area
    % Compute the area associated with each Chebyshev point using the values
    % halfway between the point below and above
    zi = zgrid_reader;
    [~, z1dc] = cheb(Nzc);
    
    % first do it on the standard interval
    % the bottom and top most pts get a half grid box
    arc(1) = 0.5*(z1dc(1)-z1dc(2));
    arc(Nzc+1) = arc(1);
    
    % over the interior pts and store the grid boxes
    arc(2:Nzc) = 0.5*(z1dc(1:end-2)-z1dc(2:end-1))+0.5*(z1dc(2:end-1)-z1dc(3:end));
    
    % now for each x point, stretch or shrink according to the local depth
    Lznow = max(zi, [], 2) - min(zi, [], 2);
    arcPhys = arc.*Lznow;
    
    % Remove values that correspond to missing z values
    arcPhys(~(zInds)) = 0;
    
    % and get a long vector of the areas and the total area
    arcPhysv = arcPhys(:);
    totar = sum(arcPhysv);
    
else
    arcPhys = ones(size(x))*(params.Lx/params.Nx)*(params.Lz/params.Nz);
    totar = sum(arcPhys(:));
end
% brutally inefficient but will work for 2D double loop
for jj = 1:Nx
    for kk = 1:Nz
        % update the corect box's total with the current area value
        myhist(var1box(jj, kk),var2box(jj, kk)) = myhist(var1box(jj, kk), var2box(jj, kk))+1*arcPhys(jj,kk);
    end
end

usp = myhist'/totar;

%% Plot up
if (nargin > 5 && isPlot) || (nargout == 0)
    isPlot = true;
else
    isPlot = false;
end
[axLab1] = get_axis_labels(var1);
[axLab2] = get_axis_labels(var2);

if isPlot
    clf;
    % Change the figure aspect ratio to taller
    %fig = gcf; %fig.Position([3 4]) = [643.2000 531.2000];
    
    %tiledlayout(4, 1, 'TileSpacing', 'compact');
    isSanityCheck = true;
    if isSanityCheck % These are the upper plots of the variables in physical space
        %figure(1)
        %ax1 = nexttile;
        ax1 = subaxis(4, 1, 1, 'MarginTop', .05); % Uncomment if using a version
        %of matlab without tiledlayouts
        pcolor(x, z, data1), shading flat;
        title(['t = ', num2str(ii)]);
        colormap(gca, cmocean('dense'));
        caxis([var1min var1max]);
        c1 = colorbar(ax1, 'Location','EastOutside');
        ylabel(c1, axLab1); ylabel('z (m)');
        axis([xlims zlims])
        xticklabels([]);
        hold on;
        plot(x(:, 1), z(:, 1), 'k-');
        
        %ax1 = nexttile;
        ax2 = subaxis(4, 1, 2, 'MarginTop', .05); % Uncomment if using a version
        %of matlab without tiledlayoutsax2 = subaxis(4, 1, 2, 'MT', 0.04);
        pcolor(x,z,data2), shading flat;
        
        if strcmpi(var2, 'vorty')
            colormap(gca, cmocean('balance'));
        else
            colormap(gca, cmocean('amp'))
        end
        ax2.Position(3) = ax1.Position(3);
        
        caxis([var2min var2max]);
        c2 = colorbar('Location','EastOutside');
        ylabel(c2, axLab2);
        xlabel('x (m)'); ylabel('z (m)');
        
        c2.Position(1) = c1.Position(1);
        ax2.Position(3) = ax1.Position(3);
        
        hold on;
        plot(x(:, 1), z(:, 1), 'k-');
        axis([xlims zlims])
    end
    
    
    ax3 = subaxis(6, 3, 1, 4, 1, 2, 'MarginBottom', 0.04, 'Holdaxis', true);
    ax4 = subaxis(6, 3, 2, 4, 1, 2, 'MarginBottom', 0.04, 'Holdaxis', true);
    ax5 = subaxis(6, 3, 2, 6, 1, 1, 'MarginBottom', 0.04, 'Holdaxis', true);
    
    % Plot the 2d QSP histogram
    axes(ax4)
    imagesc(ax4, myVar1, myVar2, log10(usp)); set(ax4, 'YDir', 'normal')
    shading flat
    yticklabels([]);
    xticklabels([]);
    colormap(gca, plasma);
    caxis([-6 -2]);
    c = colorbar; ylabel(c, 'Volume');
    axis([var1min var1max var2min var2max]);
    box on
    ax4.Position = ax4.Position;
    
    % 1D Histogram for variable 2
    axes(ax3)
    plot(ax3, sum(usp, 2)./sum(usp(:)), myVar2, 'r-');
    ylim([var2min var2max])
    xlim([0 .15]);
    ylabel(axLab2);
    ax3.Position = [ax4.Position(1)-ax5.Position(4) ax4.Position(2) ax5.Position(4) ax4.Position(4)];% plonks the axes on the edge of the QSP axis
    xticklabels([]);
    
    % 1D Histogram for variable 1
    axes(ax5)
    plot(myVar1, sum(usp)./sum(usp(:)), 'b-');
    xlim([var1min var1max])
    ylim([0 .15]);
    xlabel(axLab1);
    ax5.Position = [ax4.Position(1) ax4.Position(2)-ax5.Position(4) ax4.Position(3) ax5.Position(4)]; % plonks the axes on the edge of the QSP axis
    yticklabels([]);
    figure_print_format(gcf);
    
end
if (nargout == 0)
    clear usp
end
if (nargout > 3)
    varLims = [var2min var2max var1min var1max];
end
