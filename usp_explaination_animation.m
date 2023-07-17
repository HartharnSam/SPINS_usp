clearvars; close all; clc;
ii = 70;
var1 = 'rho'; var2 = 'KE';

phys_lims = [5.5 7 -.3 0];
var_lims = [0 1e-2 -0.0095 0.0095];
% Option to invert the selected region (show outside the rectangle)
isInvert = false;

%% Load in data
% Compute the qsp data
params = spins_params;
xlims = phys_lims([1 2]);
zlims = phys_lims([3 4]);
[qsp, myVar1, myKE, var_lims] = qsp_mapped(ii, var1, var2, phys_lims, var_lims);

close all;
var2_range = var_lims(2)-var_lims(1);
var1_range = var_lims(4)-var_lims(3);

var2_lims = [-1 -1 -1 -1 linspace(var_lims(1), var2_range/3, 30); 500 500 500 500 linspace(var_lims(2), var_lims(2)-var2_range/3, 30)];
var1_lims = [-1 -1 -1 -1 linspace(var_lims(3), var_lims(3)+var1_range/3, 30); 1 1 1 1 linspace(var_lims(4), var_lims(4)-var1_range/3, 30)];

x = xgrid_reader();
xmin_ind = nearest_index(x(:, 1), xlims(1));
xmax_ind = nearest_index(x(:, 1), xlims(2));
x = x(xmin_ind:xmax_ind, :);
z = zgrid_reader(xmin_ind:xmax_ind, []);

switch var1
    case 's'
        data1 = spins_reader_new('s',ii, xmin_ind:xmax_ind, []);
        data1 = data1.*(data1>0); % We always know for salinity -ve isn't real
    case 'rho'
        data1 = spins_reader_new('rho', ii, xmin_ind:xmax_ind, []);
    otherwise
        try
            data1 = spins_reader_new(var1, ii, xmin_ind:xmax_ind, []);
        catch
            error([var1, ' not configured']);
        end
end

switch lower(var2)
    case 'ke'
        u = spins_reader_new('u',ii, xmin_ind:xmax_ind, []);
        w = spins_reader_new('w',ii, xmin_ind:xmax_ind, []);
        data2 = 0.5*(u.^2+w.^2);
    case 'vorty'
        data2 = spins_reader_new('vorty', ii, xmin_ind:xmax_ind, []);
    case 'enstrophy'
        try
            data2 = spins_reader_new('enstrophy', ii, xmin_ind:xmax_ind, []);
        catch
            data2 = 0.5*spins_reader_new('vorty', ii, xmin_ind:xmax_ind, []).^2;
        end
    case 'diss'
        data2 = spins_reader_new('diss', ii, xmin_ind:xmax_ind, []);
        data2 = log10(data2);
    case 'speed'
        u = spins_reader_new('u',ii, xmin_ind:xmax_ind, []);
        w = spins_reader_new('w',ii, xmin_ind:xmax_ind, []);
        data2 = sqrt(u.^2 + w.^2);
    otherwise
        try
            data2 = spins_reader_new(var2, ii, xmin_ind:xmax_ind, []);
        catch
            error([var2, ' not configured']);
        end
end


vid = VideoWriter('C:\Users\samha\OneDrive - Newcastle University\02_PhD_Project\07_CanadaMixing\06_Communication\QSP2Phys.mp4', 'MPEG-4');
vid.FrameRate = 1;
vid.Quality = 100;
open(vid)


for ii = 1:length(var2_lims)


    RegOfInterest = ~((data1 >= var1_lims(1, ii)) & (data1 <= var1_lims(2, ii)) & (data2 >= var2_lims(1, ii))...
        & (data2 <= var2_lims(2, ii)));
    data1_temp = data1;
    data2_temp = data2;
    data1_temp(RegOfInterest) = NaN;
    data2_temp(RegOfInterest) = NaN;

    clf
    ax1 = subplot(3, 1, 1);
    pcolor(ax1, x, z, data1_temp); shading flat;
    colormap(gca, cmocean('dense'));
    caxis(var_lims([3 4]));
    c = colorbar('location', 'EastOutside');
    ylabel(c, var1); ylabel('z (m)');
    axis([xlims zlims])
    xticklabels([]);
    hold on;
    plot(x(:, 1), z(:, 1), 'k-');

    ax2 = subplot(3, 1, 2);
    pcolor(ax2, x, z, data2_temp); shading flat;
    if strcmpi(var2, 'vorty')
        colormap(gca, cmocean('balance'));
    else
        colormap(gca, cmocean('amp'))
    end
    caxis(var_lims([1 2]));
    c = colorbar('Location', 'EastOutside');
    ylabel(c, var2); xlabel('x (m)'); ylabel('z (m)');
    hold on;
    plot(x(:, 1), z(:, 1), 'k-');

    axis([xlims zlims])

    % Plot the QSP part
    ax3 = subplot(3, 1, 3);
    pcolor(myVar1, myKE, log10(qsp));
    shading flat;
    colormap(gca, plasma);
    caxis([-6 -2]);
    c = colorbar; ylabel(c, 'Volume');
    axis([var_lims(3) var_lims(4) var_lims(1) var_lims(2)]);
    axis square;
    box on
    hold on
    rectangle('Position', [var1_lims(1, ii) var2_lims(1, ii) (var1_lims(2, ii)-var1_lims(1,ii)) (var2_lims(2, ii)-var2_lims(1,ii))],...
        'EdgeColor', 'w');
    xlabel(var1); ylabel(var2);

    drawnow; pause(.1)
 drawnow; pause(0.1); % This seems to be needed to make sure it's done it's thing before "getting frame"
    F = getframe(gcf);
    writeVideo(vid, F);
    completion(ii, 60);
end
close(vid);
