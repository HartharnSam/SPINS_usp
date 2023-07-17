function usp_to_physical_3d(ii, var1, var2, phys_lims, var_lims, var1_lim, var2_lim)

isInvert = false;

params = spins_params;
if nargin < 4
    xlims = [params.min_x params.min_x+params.Lx];
    phys_lims = xlims;
else
    xlims = phys_lims([1 2]);
end
if nargin >=4 && numel(phys_lims) == 6
    ylims = phys_lims([3 4]);
    zlims = phys_lims([5 6]);
else
    ylims = [params.min_y params.min_y+params.Ly];
    zlims = [params.min_z params.min_z+params.Lz];
    phys_lims = [phys_lims ylims zlims];
end

if nargin>4
    [~, ~, ~, var_lims] = usp_3d(ii, var1, var2, phys_lims, var_lims);
else
    [~, ~, ~, var_lims] = usp_3d(ii, var1, var2, phys_lims);
end

%% Set the region of interest
%var1_lim = [-0.03 -0.01];%0.99999];
%var2_lim = [0 6.5]*1e-5;

%disp('Click the corners on the QSP plot to select your region of interest')
%[var1_lim, var2_lim] = ginput(2);
%var1_lim = sort(var1_lim); var2_lim = sort(var2_lim); % Sorting means the clicking can be in any 

% TODO: Add in interactive version?

%% Read in grids
x = xgrid_reader();
y = ygrid_reader();
z = zgrid_reader();

xmin_ind = nearest_index(x(:,1,  1), xlims(1));
xmax_ind = nearest_index(x(:,1, 1), xlims(2));
ymin_ind = nearest_index(y(1, :, 1), ylims(1));
ymax_ind = nearest_index(y(1, :, 1), ylims(2));
zmin_ind = nearest_index(z(1, 1, :), zlims(1));
zmax_ind = nearest_index(z(1, 1, :), zlims(2));

x = x(xmin_ind:xmax_ind, ymin_ind:ymax_ind, zmin_ind:zmax_ind);
y = y(xmin_ind:xmax_ind, ymin_ind:ymax_ind, zmin_ind:zmax_ind);
z = z(xmin_ind:xmax_ind, ymin_ind:ymax_ind, zmin_ind:zmax_ind);

slice_grids;
%[Nx, Ny, Nz] = size(x);

%% Read in data
switch var1
    case 's'
        data1 = spins_reader_new('s', ii, xmin_ind:xmax_ind, ymin_ind:ymax_ind, zmin_ind:zmax_ind);
        data1 = data1.*(data1>0);
    case 'rho'
        try 
            data1 = spins_reader_new('rho', ii, xmin_ind:xmax_ind, ymin_ind:ymax_ind, zmin_ind:zmax_ind);
        catch
            rho0 = params.rho_0;
            data1 = (eqn_of_state(spins_reader_new('t', ii, xmin_ind:xmax_ind,...
                    ymin_ind:ymax_ind, zmin_ind:zmax_ind)));
        end
    otherwise
        data1 = spins_reader_new(var1, ii, xmin_ind:xmax_ind, ymin_ind:ymax_ind, zmin_ind:zmax_ind);

end

switch lower(var2)
    case 'ke'
        u = spins_reader_new('u', ii, xmin_ind:xmax_ind, ymin_ind:ymax_ind, zmin_ind:zmax_ind);
        v = spins_reader_new('v', ii, xmin_ind:xmax_ind, ymin_ind:ymax_ind, zmin_ind:zmax_ind);
        w = spins_reader_new('w', ii, xmin_ind:xmax_ind, ymin_ind:ymax_ind, zmin_ind:zmax_ind);

        data2 = 0.5.*(u.^2 + w.^2 + v.^2);
        clear u w;
    case 'ke_v'
        v = spins_reader_new('v', ii, xmin_ind:xmax_ind, ymin_ind:ymax_ind, zmin_ind:zmax_ind);
        data2 = 0.5.*v.^2;
end

%% Calculate some statistics
data1_bar = mean(data1, 2);
data1_std = std(data1, 0, 2);

data2_bar = mean(data2, 2);
data2_std = std(data2, 0, 2);

[~, error_ind] = min(sum(sum(data1 - data1_bar, 1), 3));
data1_slice = squeeze(data1(:, error_ind, :));
data2_slice = squeeze(data2(:, error_ind, :));
%vdata_slice = squeeze(v(:, error_ind, :));

% %% Extract USP region of interest from all of the physical data
% RegOfInterest = (~((data1 >= var1_lim(1)) & (data1 <= var1_lim(2)) & (data2 >= var2_lim(1))...
%     & (data2 <= var2_lim(2))));
% if isInvert
%     data1(~RegOfInterest) = NaN;
%     data2(~RegOfInterest) = NaN;
% else
%     data1(RegOfInterest) = NaN;
%     data2(RegOfInterest) = NaN;
% end

%% Extract USP region of interest from the mean and from the slice
RegOfInterest_mean = (~((data1_bar >= var1_lim(1)) & (data1_bar <= var1_lim(2)) & (data2_bar >= var2_lim(1))...
    & (data2_bar <= var2_lim(2))));

if isInvert
    RegOfInterest = ~RegOfInterest;
end

data1_bar(RegOfInterest_mean) = NaN;
data2_bar(RegOfInterest_mean) = NaN;

%data1_std(RegOfInterest_mean) = NaN;
%data2_std(RegOfInterest_mean) = NaN;

RegOfInterest = (~((data1 >= var1_lim(1)) & (data1 <= var1_lim(2)) & (data2 >= var2_lim(1))...
    & (data2 <= var2_lim(2))));

data_trans_freq = squeeze(sum(~RegOfInterest, 2));

RegOfInterest_slice = (~((data1_slice >= var1_lim(1)) & (data1_slice <= var1_lim(2)) & (data2_slice >= var2_lim(1))...
    & (data2_slice <= var2_lim(2))));

data1_slice(RegOfInterest_slice) = NaN;
data2_slice(RegOfInterest_slice) = NaN;

%vdata_slice(RegOfInterest_slice) = NaN;

%% Plotting up
figure;
%title(['t = ', num2str(ii*params.plot_interval), 's'])
n_vert = 2;
subaxis(n_vert, 2, 1)
pcolor(xv, zv, squeeze(data1_bar)); shading flat; c = colorbar;
ylabel(c, ['$\overline{', var1, '}$'], 'interpreter', 'latex', 'FontSize', 16);
ylabel('z (m)'); 
xticklabels([]);
caxis(var_lims([3 4]));
cmocean('dense');
axis tight
title('(a)');

%subaxis(n_vert, 2, 2)
%pcolor(xv, zv, squeeze(data2_bar)); shading flat; c = colorbar;
%ylabel(c, ['$\overline{', var2, '}$'], 'interpreter', 'latex', 'FontSize', 16);
%axis tight
%yticklabels([]);xticklabels([]);
%caxis(var_lims([1 2]));
%cmocean('amp');

subaxis(n_vert, 2, 3)
pcolor(xv, zv, data_trans_freq/params.Ny *100); shading flat; c = colorbar;
ylabel(c, '% in domain');
ylabel('z (m)');
axis tight
%xticklabels([]);
caxis([0 100])
cmocean('matter');
xlabel('x (m)');
title('(c)');

%subaxis(n_vert, 2, 4)
%pcolor(xv, zv, vdata_slice.^2); shading flat; c = colorbar;
%ylabel(c, ['KE(v) (', num2str(error_ind), ')']);
%yticklabels([]);
%xticklabels([]);
%cmocean('amp');

% subaxis(n_vert, 2, 5);
% pcolor(xv, zv, data1_slice); shading flat; c = colorbar;
% ylabel(c, [var1, '(', num2str(error_ind), ')']);
% ylabel('z (m)');
% xlabel('x (m)');
% caxis(var_lims([3 4]));
% cmocean('dense');

subaxis(n_vert, 2, 2);
pcolor(xv, zv, data2_slice); shading flat; c = colorbar;
ylabel(c, ['$',var2, '(', num2str(error_ind), ')$'], 'interpreter', 'latex');
yticklabels([]);
xlabel('x (m)');
caxis(var_lims([1 2]));
cmocean('amp');
title('(b)')

%subaxis(n_vert, 2, 6);
%pcolor(xv, zv, data2_slice); shading flat; c = colorbar; 
%ylabel(c, [var2, '(', num2str(error_ind), ')']);
%yticklabels([]);
%xlabel('x (m)');
%caxis(var_lims([1 2]));
%cmocean('amp');

%subaxis(5, 2, 1, 4, 2, 2)
%pcolor(myVar1, myKE, log10(qsp)); shading flat; colorbar;
%axis square
%ylabel(var2); xlabel(var1);
%colormap(gca, plasma);
%hold on
%rectangle('Position', [var1_lim(1) var2_lim(1) diff(var1_lim) diff(var2_lim)], 'EdgeColor', 'w');
