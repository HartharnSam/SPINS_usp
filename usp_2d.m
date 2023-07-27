function [usp, myVar1, myVar2, var_lims] = usp_2d(ti, var1, var2, spat_lims, var_lims, doPlot)
%USP_2D - produces USP, or paired histograms to tell us where two
%variables overlap, so we can find out if, when and where we get
%combinations of two variables.
%
% Syntax:  [usp, myVar1, myVar2, var_lims] = usp_2d(ti, var1, var2, spat_lims, var_lims, doPlot
%
% Inputs:
%    ti - Simulation timestep to output for
%    var1 - variable to compare to var2 (on x axis)
%    var2 - Variable to compare to var1 (on y axis)
%    spat_lims - [optional] Spatial Region of physical space [xmin xmax zmin zmax]
%    specify the x limits. Defaults to full size of tank
%    var_lims - [optional] limits of variables to investigate as [var2min var2max var1min var2max]
%
%   [NOTE] - var lims are in the order var2, var1 - it seems
%   counterintuative, but more likely to set only the y limits than only x!
%
% Outputs:
%    usp - Histogram of combined volume
%    myVar1 - variable values for each histogram box
%    myVar2 - ""
%    var_lims - Limits of the second variable calculated by the script
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
%
% Dev notes:
% - I think this won't work for anything that isn't FOURIER in x, cheb in
% z. Relatively minor changes needed if they are simply the other way
% around, more if they are both fourier?
%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

%% read in the grids 
% & cut down to size
params = spins_params;
if nargin < 4 || isempty(spat_lims)
    xlims = [params.min_x params.min_x+params.Lx];
else
    xlims = spat_lims([1 2]);
end
if nargin >=4 && numel(spat_lims) == 4
    zlims = spat_lims([3 4]);
else
    zlims = [params.min_z params.min_z+params.Lz];
end

x = xgrid_reader();
z = zgrid_reader();

xmin_ind = nearest_index(x(:, 1), xlims(1));
xmax_ind = nearest_index(x(:, 1), xlims(2));
if strcmpi(params.mapped_grid, 'true')
    z_inds = [];
    x = x(xmin_ind:xmax_ind, :);
    z = z(xmin_ind:xmax_ind, :);


else
    zmin_ind = nearest_index(z(1, :), zlims(1));
    zmax_ind = nearest_index(z(1, :), zlims(2));
    z_inds = zmin_ind:zmax_ind;
    x = x(xmin_ind:xmax_ind, z_inds);
    z = z(xmin_ind:xmax_ind, z_inds);
end
% get the size of the grids
[Nx, Nz] = size(x);

% parameter for the chebyshev grid which are on [-1,1]
Nzc = Nz-1;
%% read in data
switch lower(var1)
    case 's'
        data1 = spins_reader_new('s', ti, xmin_ind:xmax_ind, z_inds);
        data1 = data1.*(data1>0);
    case 'enstrophy'
        try
            data1 = spins_reader_new('enst', ti, xmin_ind:xmax_ind, z_inds);
        catch
            data1 = 0.5*spins_reader_new('vorty', ti, xmin_ind:xmax_ind, z_inds).^2;
        end
    case 'rho'
        try
            data1 = spins_reader_new('rho', ti, xmin_ind:xmax_ind, z_inds);
        catch
            rho0 = params.rho_0;
            data1 = (eqn_of_state(spins_reader_new('t', ti, xmin_ind:xmax_ind,...
                z_inds), 0));
        end
    otherwise
        try
            data1 = spins_reader_new(var1, ti, xmin_ind:xmax_ind, z_inds);
        catch
            error([var1, ' not configured']);
        end
end

switch lower(var2)
    case 'ke'
        u = spins_reader_new('u', ti, xmin_ind:xmax_ind, z_inds);
        w = spins_reader_new('w', ti, xmin_ind:xmax_ind, z_inds);
        data2 = 0.5*(u.^2 + w.^2);
        clear u w
    case 'vorty'
        data2 = spins_reader_new('vorty', ti, xmin_ind:xmax_ind, z_inds);
    case 'enstrophy'
        try
            data2 = spins_reader_new('enst', ti, xmin_ind:xmax_ind, z_inds);
        catch
            data2 = 0.5*spins_reader_new('vorty', ti, xmin_ind:xmax_ind, z_inds).^2;
        end
    case 'diss'
        data2 = spins_reader_new('diss', ti, xmin_ind:xmax_ind, z_inds);
        data2 = log10(data2);
    case 'speed'
        u = spins_reader_new('u',ti, xmin_ind:xmax_ind, z_inds);
        w = spins_reader_new('w',ti, xmin_ind:xmax_ind, z_inds);
        data2 = sqrt(u.^2 + w.^2);
        clear u w
    otherwise
        try
            data2 = spins_reader_new(var2, ti, xmin_ind:xmax_ind, z_inds);
        catch
            error([var2, ' not configured']);
        end
end

%% Sort data into the histogram "boxes"
numpts = 50;
% for variable on x
if nargin > 4 && numel(var_lims) == 4
    var1min = var_lims(3); var1max = var_lims(4);
else % Use defaults for the var1 limits
    if strcmpi(var1, 'rho')
        var1min = -params.delta_rho/2;
        var1max = -var1min;
    else
        var1min = min(data1(:));
        var1max = max(data1(:));
    end
end
% Cap the values at the limits - ignoring these values by NaN'ing them
% would be a reasonable choice, but requires some more coding to make it
% work
data1(data1<var1min) = var1min;
data1(data1>var1max) = var1max;

ranges = var1max-var1min;
dVar1 = ranges/(numpts-1);
myVar1 = var1min+(0.5:numpts-0.5)*dVar1;

% For variable on y (KE)
if nargin <= 4 || isempty(var_lims) % Use defaults
    var2min = min(data2(:));
    var2max = max(data2(:));
else
    var2min = var_lims(1);
    var2max = var_lims(2);
end
data2(data2<var2min) = var2min;
data2(data2>var2max) = var2max;

rangeVar2 = var2max-var2min;
dVar2 = rangeVar2/(numpts-1);
myVar2 = var2min+(0.5:numpts-0.5)'*dVar2;

myhist = zeros(numpts,numpts);

%figure out which box coordinate you are in
var1box = ceil((data1-var1min)/dVar1);
var1box = var1box+1*(data1==var1min);

var2box = ceil((data2-var2min)/dVar2);
var2box = var2box+1*(data2==var2min);

if strcmpi(params.mapped_grid, 'true')
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
    arcphys = arc.*Lznow;

    % Remove values that correspond to missing z values
    arcphys(~(z_inds)) = 0;

    % and get a long vector of the areas and the total area
    arcphysv = arcphys(:);
    totar = sum(arcphysv);

else
    arcphys = ones(size(x))*(params.Lx/params.Nx)*(params.Lz/params.Nz);
    totar = sum(arcphys(:));
end
% brutally inefficient but will work for 2D double loop
for i = 1:Nx
    for jj = 1:Nz
        % update the corect box's total with the current area value
        myhist(var1box(i, jj),var2box(i, jj)) = myhist(var1box(i, jj), var2box(i, jj))+1*arcphys(i,jj);
    end
end
usp = myhist'/totar;

%% Plot up
if nargin > 5 && doPlot
    isPlot = true;
else
    isPlot = false;
end
if isPlot
    clf;
    % Change the figure aspect ratio to taller
    fig = gcf; %fig.Position([3 4]) = [643.2000 531.2000];

    %tiledlayout(4, 1, 'TileSpacing', 'compact');
    isSanityCheck = true;
    if isSanityCheck % These are the upper plots of the variables in physical space
        figure(1)
        %ax1 = nexttile;
        ax1 = subaxis(4, 1, 1, 'MT', .05); % Uncomment if using a version
        %of matlab without tiledlayouts
        pcolor(x, z, data1), shading flat;
        title(['t = ', num2str(ti)]);
        colormap(gca, cmocean('dense'));
        caxis([var1min var1max]);
        c1 = colorbar(ax1, 'Location','EastOutside');
        ylabel(c1, var1); ylabel('z (m)');
        axis([xlims zlims])
        xticklabels([]);
        hold on;
        plot(x(:, 1), z(:, 1), 'k-');

        %ax1 = nexttile;
        ax2 = subaxis(4, 1, 2, 'MT', .05); % Uncomment if using a version
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
        ylabel(c2, var2);
         xlabel('x (m)'); ylabel('z (m)');

        c2.Position(1) = c1.Position(1);
        ax2.Position(3) = ax1.Position(3);

        hold on;
        plot(x(:, 1), z(:, 1), 'k-');
        axis([xlims zlims])
    end


    ax3 = subaxis(6, 3, 1, 4, 1, 2, 'MB', 0.04, 'Holdaxis');
    ax4 = subaxis(6, 3, 2, 4, 1, 2, 'MB', 0.04, 'Holdaxis');
    ax5 = subaxis(6, 3, 2, 6, 1, 1, 'MB', 0.04, 'Holdaxis');

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
    ylabel(var2);
    ax3.Position = [ax4.Position(1)-ax5.Position(4) ax4.Position(2) ax5.Position(4) ax4.Position(4)];% plonks the axes on the edge of the QSP axis
    xticklabels([]);

    % 1D Histogram for variable 1
    axes(ax5)
    plot(myVar1, sum(usp)./sum(usp(:)), 'b-');
    xlim([var1min var1max])
    ylim([0 .15]);
    xlabel(var1);
    ax5.Position = [ax4.Position(1) ax4.Position(2)-ax5.Position(4) ax4.Position(3) ax5.Position(4)]; % plonks the axes on the edge of the QSP axis
    yticklabels([]);
    figure_print_format(gcf);

end
if nargout > 3
    var_lims = [var2min var2max var1min var1max];
end