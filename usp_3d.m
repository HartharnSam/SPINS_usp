function [usp, myVar1, myVar2, varLims] = usp_3d(ii, var1, var2, spatLims, varLims, opts)
%USP_3D - produces USP, or paired histograms to tell us where two
%variables overlap, so we can find out if, when and where we get
%combinations of two variables. For 3D unmapped simulations
%
% Syntax:  [usp, myVar1, myVar2, varLims] = usp_3d(ii, var1, var2, physLims, varLims)
%
% Inputs:
%    ii - Simulation timestep to output for
%    var1 - variable to compare to var2 (on x axis)
%    var2 - Variable to compare to var1 (on y axis)
%    physLims - [optional] Region of physical space [xmin xmax ymin ymax zmin zmax], can only
%    specify the x limits. Defaults to full size of tank
%    varLims - [optional] limits of variables to investigate as [var2min var2max var1min var2max]
%
%   [NOTE] - var lims are in the order var2, var1 - it seems
%   counterintuative, but more likely to set only the y limits than only x!
%
% Outputs:
%    usp - Histogram of combined volume
%    myVar1 - variable values for each histogram box
%    myVar2 - ""
%    varLims - Limits of the second variable calculated by the script
%
% Other m-files required: xgrid_reader, zgrid_reader, spins_params,
% spins_reader_new, nearest_index, cheb, figure_print_format,
% plasma, subaxis
% Subfunctions: none
% MAT-files required: none
%
% See also: usp_to_physical, spins_QSP_csv
% Author:
% Department of Applied Mathematics, University of Waterloo
% email address:
% GitHub:
% 14-Jun-2022; Last revision: 16-Aug-2024
% MATLAB Version: 9.12.0.1956245 (R2022a) Update 2
arguments
    ii (1, 1) uint16
    var1 (1, 1) string
    var2 (1, 1) string
    spatLims (:, 1) double = []
    varLims (:, 1) double = []
    opts.data1 (:, :) double = []
    opts.data2 (:, :) double = []
end
clf
%% Read in the grids
% & cut down size as required
params = spins_params;
if isempty(spatLims)
    xlims = [0 params.Lx]+params.min_x;
else
    xlims = spatLims([1 2]);
end
if (numel(spatLims) == 6)
    ylims = spatLims([3 4]);
    zlims = spatLims([5 6]);
else
    ylims = [0 params.Ly]+params.min_y;
    zlims = [0 params.Lz]+params.min_z ;
end

x = xgrid_reader();
y = ygrid_reader();
z = zgrid_reader();

xminInd = nearest_index(x(:,1,  1), xlims(1));
xmaxInd = nearest_index(x(:,1, 1), xlims(2));
yminInd = nearest_index(y(1, :, 1), ylims(1));
ymaxInd = nearest_index(y(1, :, 1), ylims(2));
zminInd = nearest_index(z(1, 1, :), zlims(1));
zmaxInd = nearest_index(z(1, 1, :), zlims(2));

x = x(xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
y = y(xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
z = z(xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);

[xv, zv] = slice_grids(x, y, z);
[Nx, Ny, Nz] = size(x);

%% Read in data
if isempty(opts.data1)
    data1 = get_spins_3Ddata(var1, ii, xminInd, xmaxInd, yminInd, ymaxInd, zInds);
else
    data1 = opts.data1;
end
if isempty(opts.data2)
    data2 = get_spins_3Ddata(var2, ii, xminInd, xmaxInd, yminInd, ymaxInd, zInds);
else
    data2 = opts.data2;
end
%% Sort the data into the histogram boxes
numpts = 50;
if (nargin > 4) && (numel(varLims) == 4)
    var1min = varLims(3); var1max = varLims(4);
else % Use default var1 limits
    var1min = min(data1(:));
    var1max = max(data1(:));
end
data1(data1 < var1min) = var1min;
data1(data1 > var1max) = var1max;

if nargin > 4
    var2min = varLims(1);
    var2max = varLims(2);
else
    var2min = min(data2(:));
    var2max = max(data2(:));
end
data2(data2 < var2min) = var2min;
data2(data2 > var2max) = var2max;

rangeVar1 = var1max-var1min;
rangeVar2 = var2max-var2min;

dVar1 = rangeVar1/(numpts-1);
dVar2 = rangeVar2/(numpts-1);

var1box = ceil((data1-var1min)/dVar1);
var1box = var1box+1*(data1 == var1min);

var2box = ceil((data2-var2min)/dVar2);
var2box = var2box+1*(data2 == var2min);

myVar1 = var1min+(0.5:numpts-0.5)'*dVar1;
myVar2 = var2min+(0.5:numpts-0.5)'*dVar2;

myhist = zeros(numpts, numpts);

isCheb = isequal(params.mapped_grid, 'true') || isequal(params.type_z, 'NO_SLIP');
assert(~isCheb, "Currently not set up for chebyshev grids, see usp_2d for how to implement");

for jj = 1:Nx
    for kk = 1:Ny
        for ll = 1:Nz
            myhist(var1box(jj, kk, ll), var2box(jj, kk, ll)) = myhist(var1box(jj, kk, ll), var2box(jj, kk, ll))+1;
        end
    end
end
usp = myhist'; %TODO; weight by total area

%% Plot the data
isPlot = true; % speeds up processing with no graphical outputs needed
if isPlot
    clf;
    figure(1)
    
    subaxis(4, 1, 1, 'MarginTop', 0.04);
    pcolor(xv, zv, squeeze(data1(:, 1, :))); shading flat;
    title(['t = ', num2str(ii*params.plot_interval)]);
    colormap(gca, cmocean('dense'));
    caxis([var1min var1max]);
    c = colorbar('Location', 'EastOutside');
    ylabel(c, var1); ylabel('z (m)');
    axis([xlims zlims])
    
    subaxis(4, 1, 2, 'MarginTop', 0.04);
    pcolor(xv, zv, squeeze(data2(:, 1, :))); shading flat;
    colormap(gca, cmocean('amp'))
    caxis([var2min var2max]);
    c = colorbar('Location', 'EastOutside');
    ylabel(c, var2); ylabel('z (m)');
    axis([xlims zlims])
    
    ax3 = subaxis(6, 3, 1, 4, 1, 2, 'MB', 0.04);
    ax4 = subaxis(6, 3, 2, 4, 1, 2, 'MB', 0.04);
    ax5 = subaxis(6, 3, 2, 6, 1, 1, 'MB', 0.04);
    
    % plot 2D QSP histogram
    axes(ax4)
    imagesc(ax4, myVar1, myVar2, log10(usp)); shading flat; set(gca, 'YDir', 'normal')
    yticklabels([]); xticklabels([]);
    colormap(gca, plasma);
    %caxis([-6 -2]);
    c = colorbar; ylabel(c, 'Volume')
    axis([var1min var1max var2min var2max]);
    box on
    ax4.Position = ax4.Position;
    
    axes(ax3)
    plot(ax3, sum(usp, 2)./sum(usp(:)), myVar2, 'r-');
    ylim([var2min var2max])
    xlim([0 .15]);
    ylabel(var2);
    ax3.Position = [ax4.Position(1)-ax5.Position(4) ax4.Position(2) ax5.Position(4) ax4.Position(4)];% plonks the axes on the edge of the USP axis
    xticklabels([]);
    
    % 1D Histogram for variable 1
    axes(ax5)
    plot(myVar1, sum(usp)./sum(usp(:)), 'b-');
    xlim([var1min var1max])
    ylim([0 .15]);
    xlabel(var1);
    ax5.Position = [ax4.Position(1) ax4.Position(2)-ax5.Position(4) ax4.Position(3) ax5.Position(4)]; % plonks the axes on the edge of the USP axis
    yticklabels([]);
end
if (nargout == 0)
    clear usp
end
if (nargout > 3)
    varLims = [var2min var2max var1min var1max];
end


