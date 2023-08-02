function usp_to_physical_3d(ii, var1, var2, physLims, varLims, var1Lim, var2Lim)
%USP_TO_PHYSICAL_3D - Relate region of usp graph to physical space for 3D sims.
% Can be used interactively by clicking a region, or by setting it in the
% code
%
% Syntax:  usp_to_physical_3d(ii, var1, var2, spatLims, varLims, region)
%
% Inputs:
%    ii - Simulation timestep to output for
%    var1 - variable to compare to var2 (on x axis)
%    var2 - Variable to compare to var1 (on y axis)
%    spatLims - [optional] Spatial region of physical space [xmin xmax ymin ymax zmin zmax]
%       optionally, only set the x limits. Defaults to full size of tank
%    varLims - [optional] realistic limits of the variables to investigate as [var2min var2max var1min var2max]
%    region - [optional] USP region of interest to display data for
%
% Other m-files required: usp_3d, spins_params, xgrid_reader,
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
isInvert = false;

params = spins_params;
if nargin < 4
    xlims = [params.min_x params.min_x+params.Lx];
    physLims = xlims;
else
    xlims = physLims([1 2]);
end
if nargin >= 4 && numel(physLims) == 6
    ylims = physLims([3 4]);
    zlims = physLims([5 6]);
else
    ylims = [params.min_y params.min_y+params.Ly];
    zlims = [params.min_z params.min_z+params.Lz];
    physLims = [physLims ylims zlims];
end

if nargin > 4
    [~, ~, ~, varLims] = usp_3d(ii, var1, var2, physLims, varLims);
else
    [~, ~, ~, varLims] = usp_3d(ii, var1, var2, physLims);
end

%% Set the region of interest
%var1Lim = [-0.03 -0.01];%0.99999];
%var2Lim = [0 6.5]*1e-5;

%disp('Click the corners on the QSP plot to select your region of interest')
%[var1Lim, var2Lim] = ginput(2);
%var1Lim = sort(var1Lim); var2Lim = sort(var2Lim); % Sorting means the clicking can be in any

% TODO: Add in interactive version?

%% Read in grids
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

slice_grids;
%[Nx, Ny, Nz] = size(x);

%% Read in data
switch var1
    case 's'
        data1 = spins_reader_new('s', ii, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        data1 = data1.*(data1 > 0);
    case 'rho'
        try
            data1 = spins_reader_new('rho', ii, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        catch
            %rho0 = params.rho_0;
            data1 = (eqn_of_state(spins_reader_new('t', ii, xminInd:xmaxInd,...
                yminInd:ymaxInd, zminInd:zmaxInd)));
        end
    otherwise
        data1 = spins_reader_new(var1, ii, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        
end

switch lower(var2)
    case 'ke'
        u = spins_reader_new('u', ii, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        v = spins_reader_new('v', ii, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        w = spins_reader_new('w', ii, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        
        data2 = 0.5.*(u.^2 + w.^2 + v.^2);
        clear u w;
    case 'ke_v'
        v = spins_reader_new('v', ii, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        data2 = 0.5.*v.^2;
end

%% Calculate some statistics
data1Bar = mean(data1, 2);
data1STDEV = std(data1, 0, 2);

data2Bar = mean(data2, 2);
data2STDEV = std(data2, 0, 2);

[~, errorInd] = min(sum(sum(data1 - data1Bar, 1), 3));
data1Slice = squeeze(data1(:, errorInd, :));
data2Slice = squeeze(data2(:, errorInd, :));
%vdata_slice = squeeze(v(:, errorInd, :));

% %% Extract USP region of interest from all of the physical data
% regionOfInterest = (~((data1 >= var1Lim(1)) & (data1 <= var1Lim(2)) & (data2 >= var2Lim(1))...
%     & (data2 <= var2Lim(2))));
% if isInvert
%     data1(~regionOfInterest) = NaN;
%     data2(~regionOfInterest) = NaN;
% else
%     data1(regionOfInterest) = NaN;
%     data2(regionOfInterest) = NaN;
% end

%% Extract USP region of interest from the mean and from the slice
regionOfInterestMean = (~((data1Bar >= var1Lim(1)) & (data1Bar <= var1Lim(2)) & (data2Bar >= var2Lim(1))...
    & (data2Bar <= var2Lim(2))));

if isInvert
    regionOfInterest = ~regionOfInterest;
end

data1Bar(regionOfInterestMean) = NaN;
data2Bar(regionOfInterestMean) = NaN;

%data1STDEV(regionOfInterestMean) = NaN;
%data2STDEV(regionOfInterestMean) = NaN;

regionOfInterest = (~((data1 >= var1Lim(1)) & (data1 <= var1Lim(2)) & (data2 >= var2Lim(1))...
    & (data2 <= var2Lim(2))));

dataTransFreq = squeeze(sum(~regionOfInterest, 2));

regionOfInterest_slice = (~((data1Slice >= var1Lim(1)) & (data1Slice <= var1Lim(2)) & (data2Slice >= var2Lim(1))...
    & (data2Slice <= var2Lim(2))));

data1Slice(regionOfInterest_slice) = NaN;
data2Slice(regionOfInterest_slice) = NaN;

%vdata_slice(regionOfInterest_slice) = NaN;

%% Plotting up
figure;
%title(['t = ', num2str(ii*params.plot_interval), 's'])
nVert = 2;
subaxis(nVert, 2, 1)
pcolor(xv, zv, squeeze(data1Bar)); shading flat; c = colorbar;
ylabel(c, ['$\overline{', var1, '}$'], 'interpreter', 'latex', 'FontSize', 16);
ylabel('z (m)');
xticklabels([]);
caxis(varLims([3 4]));
cmocean('dense');
axis tight
title('(a)');

%subaxis(nVert, 2, 2)
%pcolor(xv, zv, squeeze(data2Bar)); shading flat; c = colorbar;
%ylabel(c, ['$\overline{', var2, '}$'], 'interpreter', 'latex', 'FontSize', 16);
%axis tight
%yticklabels([]);xticklabels([]);
%caxis(varLims([1 2]));
%cmocean('amp');

subaxis(nVert, 2, 3)
pcolor(xv, zv, dataTransFreq/params.Ny *100); shading flat; c = colorbar;
ylabel(c, '% in domain');
ylabel('z (m)');
axis tight
%xticklabels([]);
caxis([0 100])
cmocean('matter');
xlabel('x (m)');
title('(c)');

%subaxis(nVert, 2, 4)
%pcolor(xv, zv, vdata_slice.^2); shading flat; c = colorbar;
%ylabel(c, ['KE(v) (', num2str(errorInd), ')']);
%yticklabels([]);
%xticklabels([]);
%cmocean('amp');

% subaxis(nVert, 2, 5);
% pcolor(xv, zv, data1Slice); shading flat; c = colorbar;
% ylabel(c, [var1, '(', num2str(errorInd), ')']);
% ylabel('z (m)');
% xlabel('x (m)');
% caxis(varLims([3 4]));
% cmocean('dense');

subaxis(nVert, 2, 2);
pcolor(xv, zv, data2Slice); shading flat; c = colorbar;
ylabel(c, ['$',var2, '(', num2str(errorInd), ')$'], 'interpreter', 'latex');
yticklabels([]);
xlabel('x (m)');
caxis(varLims([1 2]));
cmocean('amp');
title('(b)')

%subaxis(nVert, 2, 6);
%pcolor(xv, zv, data2Slice); shading flat; c = colorbar;
%ylabel(c, [var2, '(', num2str(errorInd), ')']);
%yticklabels([]);
%xlabel('x (m)');
%caxis(varLims([1 2]));
%cmocean('amp');

%subaxis(5, 2, 1, 4, 2, 2)
%pcolor(myVar1, myKE, log10(qsp)); shading flat; colorbar;
%axis square
%ylabel(var2); xlabel(var1);
%colormap(gca, plasma);
%hold on
%rectangle('Position', [var1Lim(1) var2Lim(1) diff(var1Lim) diff(var2Lim)], 'EdgeColor', 'w');
