function [usp, myVar1, myVar2, varLims] = usp_3d(ti, var1, var2, physLims, varLims)
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
% 14-Jun-2022; Last revision: 20-Jul-2023
% MATLAB Version: 9.12.0.1956245 (R2022a) Update 2

clf
%% Read in the grids
% & cut down size as required
params = spins_params;
if nargin < 4 || isempty(physLims)
    xlims = [params.min_x params.min_x+params.Lx];
else
    xlims = physLims([1 2]);
end
if nargin >= 4 && numel(physLims) == 6
    ylims = physLims([3 4]);
    zlims = physLims([5 6]);
else
    ylims = [params.min_y params.min_y+params.Ly];
    zlims = [params.min_z params.min_z+params.Lz];
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
switch var1
    case 's'
        data1 = spins_reader_new('s', ti, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        data1 = data1.*(data1 > 0);
    case 'rho'
        try
            data1 = spins_reader_new('rho', ti, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        catch
            %rho0 = params.rho_0;
            data1 = (eqn_of_state(spins_reader_new('t', ti, xminInd:xmaxInd,...
                yminInd:ymaxInd, zminInd:zmaxInd)));
        end
    case 'enstrophy'
        try
            data1 = spins_reader_new('enst', ti, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        catch
            data1 = 0.5*spins_reader_new('vorty', ti, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd).^2;
        end
    otherwise
        data1 = spins_reader_new(var1, ti, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        
end

switch lower(var2)
    case 'ke'
        u = spins_reader_new('u', ti, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        v = spins_reader_new('v', ti, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        w = spins_reader_new('w', ti, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        
        data2 = 0.5.*(u.^2 + w.^2 + v.^2);
        clear u v w;
    case 'ke_v'
        v = spins_reader_new('v', ti, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        data2 = 0.5.*(v.^2);
        clear v;
    case 'speed'
        u = spins_reader_new('u', ti, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        v = spins_reader_new('v', ti, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        w = spins_reader_new('w', ti, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        
        data2 = sqrt(u.^2 + w.^2 + v.^2);
        clear u v w;
    case 'diss'
        data2 = log10(spins_reader_new('diss', ti, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd));
        
    otherwise
        data2 = spins_reader_new(var2, ti, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        % TODO: Add in vorticities, enstrophy?
end

%% Sort the data into the histogram boxes
numpts = 50;
if nargin > 4 && numel(varLims) == 4
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

for ii = 1:Nx
    for jj = 1:Ny
        for kk = 1:Nz
            myhist(var1box(ii, jj, kk), var2box(ii, jj, kk)) = myhist(var1box(ii, jj, kk), var2box(ii, jj, kk))+1;
        end
    end
end
usp = myhist'; %TODO; weight by total area

%% Plot the data
isPlot = true; % speeds up processing with no graphical outputs needed
if isPlot
    clf;
    figure(1)
    
    subaxis(4, 1, 1, 'MT', 0.04);
    pcolor(xv, zv, squeeze(data1(:, 1, :))); shading flat;
    title(['t = ', num2str(ti*params.plot_interval)]);
    colormap(gca, cmocean('dense'));
    caxis([var1min var1max]);
    c = colorbar('Location', 'EastOutside');
    ylabel(c, var1); ylabel('z (m)');
    axis([xlims zlims])
    
    subaxis(4, 1, 2, 'MT', 0.04);
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
if nargout > 3
    varLims = [var2min var2max var1min var1max];
end


