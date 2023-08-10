function usp_phys_plotter_comb(data_directory, times, var_lims, var_ROI, xlims, ylims, output_fnm, isFinal, sptitle)
%USP_PHYS_PLOTTER Summary of this function goes here
%   Detailed explanation goes here
orig_dir = cd;
addpath(genpath('../03_usp_matlab/'));

tlo = tiledlayout(3, length(data_directory),  'TileSpacing', 'compact', 'Padding', 'Compact', 'TileIndexing', 'ColumnMajor');
for ii = 1:length(data_directory)
    cd(data_directory{ii});
    [x, z] = spinsgrid2d;
    params = spins_params;

    [usp, myvar1, myvar2] = usp_2d(times(ii, end), 'rho', 'enstrophy', xlims(ii, :), var_lims(ii, :), false);
    ROI = usp_to_physical(times(ii, end), 'rho', 'enstrophy', [xlims(ii, :) ylims(ii, :)], var_lims(ii,:), var_ROI(ii, :));
    rho = spins_reader_new('rho', times(ii, end));
    %enst = (spins_reader_new('vorty', times(ii, end)).^2)*.5;
    
    nexttile;
    %pcolor(x, z, (rho*params.rho_0)+params.rho_0);cmocean('dense')
    pcolor(x, z, (rho*params.rho_0));cmocean('dense')
    %var_lims(ii, [3 4]) = round((var_lims(ii, [3 4])*params.rho_0)+params.rho_0);
    var_lims(ii, [3 4]) = round((var_lims(ii, [3 4])*params.rho_0));
    %var_ROI(ii, [1 2]) = var_ROI(ii, [1 2])*params.rho_0 + params.rho_0;
    var_ROI(ii, [1 2]) = var_ROI(ii, [1 2])*params.rho_0;

    colormap(gca, cmocean('dense'));
    caxis(var_lims(ii, [3 4]));
    if ii == length(data_directory)
        c = colorbar('location', 'EastOutside');
        ylabel(c, "$\rho' (kg m^{-3})$"); yticklabels([]);
    elseif ii == 1
        ylabel('z (m)');
    else
        yticklabels([]);
    end
    axis([xlims(ii, :) ylims(ii, :)]);
    xticklabels([]);
    hold on
    plot(x(:, 1), z(:, 1), 'k-');
    title(sptitle{ii})
    ntitle(subplot_labels(ii, 2), 'location', 'northwest', 'cover', true)

    nexttile;
    pcolor(ROI.x, ROI.z, double(~ROI.region));
    caxis([0 1]); cmocean('amp');
    if ii == 1
        ylabel('$z (m)$');
    else
        yticklabels([]);
    end
    xlabel('$x (m)$')
    hold on
    plot(x(:, 1), z(:, 1), 'k-');
    axis([xlims(ii, :) ylims(ii, :)]);
    ntitle(subplot_labels(ii+length(data_directory), 2), 'location', 'northwest', 'cover', true)

    nexttile;
    %imagesc((myvar1*params.rho_0)+params.rho_0, myvar2, log10(usp)); set(gca, 'YDir', 'normal');
    imagesc((myvar1*params.rho_0), myvar2, log10(usp)); set(gca, 'YDir', 'normal');

    colormap(gca, 'plasma'); caxis([-6 -2]);
    if ii == length(data_directory)
        c = colorbar('location', 'EastOutside');
        ylabel(c, '$\log_{10}(W)$');
    elseif ii == 1
        ylabel('$\Omega (s^{-2})$');
    end

    axis(var_lims(ii, [3 4 1 2]));
    hold on
    rectangle('Position', [var_ROI(ii, 1) var_ROI(ii, 3) var_ROI(ii, 2)-var_ROI(ii, 1) var_ROI(ii, 4)-var_ROI(ii, 3)], 'EdgeColor', 'w', 'LineWidth', 1, 'LineStyle', '-.');
    box on; axis square;
    xticks(round([var_lims(ii, 3), mean(var_lims(ii, [3 4])), var_lims(ii, 4)]))
    xlabel("$\rho' (kg m^{-3})$")
    ntitle(subplot_labels(ii+length(data_directory)*2, 2), 'location', 'north', 'color', 'w')

    cd(orig_dir);
end
if isFinal
    fig = gcf;
    figure_print_format(fig);

    set(gcf, 'PaperUnits', 'centimeters', 'Units', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 27.2 18.5], 'Position', [0 0 27.2 18.5]);
    % IF YOUR MATLAB VERSION IS BEFORE R2020a COMMENT THESE LINES
    exportgraphics(fig, ['../02_Raster/', output_fnm, '.png'], 'Resolution', 300);
    exportgraphics(fig, ['../03_Vector/', output_fnm, '.eps']);
end