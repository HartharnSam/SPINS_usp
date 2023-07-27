function usp_phys_plotter(data_directory, times, var_lims, var_ROI, xlims, ylims, output_fnm, plots, isFinal)
%USP_PHYS_PLOTTER Summary of this function goes here
%   Detailed explanation goes here
orig_dir = cd;
addpath(genpath('../03_usp_matlab/'));
cd(data_directory);
tiledlayout(3, length(times), 'TileSpacing', 'compact', 'Padding', 'Compact', 'TileIndexing', 'columnmajor')
[x, z] = spinsgrid2d;
params = spins_params;
for tt = 1:length(times)
    %% Load and plot Density
    rho = spins_reader_new('rho', times(tt));
    if plots == 1
        nexttile;
        %pcolor(x, z, (rho*params.rho_0)+params.rho_0);
        pcolor(x, z, (rho*params.rho_0));
        xlim(xlims); ylim(ylims);

        %caxis(gca, (var_lims([3,4])*params.rho_0)+params.rho_0);
        caxis(gca, (var_lims([3,4])*params.rho_0));
        colormap(gca, cmocean('dense'));
        if tt == length(times)
            c = colorbar('Location','EastOutside');
            ylabel(c, "$\rho' (kg m^{-3})$");
            yticklabels([]);
        elseif tt == 1
            ylabel('$z (m)$')
        else
            yticklabels([]);

        end
        xticklabels([]);
        title(subplot_labels(tt, 2))
    end
    %% Load and plot enstrophy
    enst = (spins_reader_new('vorty', times(tt)).^2)*.5;
    %ke = 0.5*(spins_reader_new('u', times(tt)).^2 + spins_reader_new('w', times(tt)).^2);
    if plots == 1
        nexttile;
        pcolor(x, z, enst);
        xlim(xlims);ylim(ylims);
        caxis(gca, var_lims([1,2]));
        colormap(gca, cmocean('amp'));
        if tt == length(times)
            c = colorbar('Location', 'EastOutside');
            ylabel(c, '$\Omega (s^{-2})$');
            yticklabels([]);
        elseif tt == 1
            ylabel('$z (m)$')
        else
            yticklabels([]);

        end
        xlabel('$x (m)$');
        title(subplot_labels(tt+length(times), 2));
    end
    %% Calculate & Plot USP
    [usp, myvar1, myvar2] = usp_2d(times(tt), 'rho', 'enstrophy', xlims, var_lims);
    if plots == 1
        nexttile;
        %tmp_var_lims([3 4]) = round((var_lims([3 4])*params.rho_0)+params.rho_0);
        tmp_var_lims([3 4]) = round((var_lims([3 4])*params.rho_0));

        tmp_var_lims([1 2]) = var_lims([1 2]);
        %imagesc((myvar1*params.rho_0)+params.rho_0, myvar2, log10(usp)); set(gca, 'YDir', 'normal');
        imagesc((myvar1*params.rho_0), myvar2, log10(usp)); set(gca, 'YDir', 'normal');

        colormap(gca, plasma);
        caxis(gca, [-6 -2]);
        if tt == length (times)
            c = colorbar; ylabel(c, '$\log_{10}(W)$');
            yticklabels([]);
        elseif tt == 1
            ylabel('$\Omega (s^{-2})$');
        else
            yticklabels([]);
        end
        xlabel("$\rho' (kg m^{-3})$")
        axis(gca, [tmp_var_lims(3) tmp_var_lims(4) var_lims(1) var_lims(2)]);
        xticks([tmp_var_lims(3) (tmp_var_lims(3)+tmp_var_lims(4))/2 tmp_var_lims(4)]);
        box on
        axis square
        if tt == 1
            tmpVarROI = var_ROI; tmpVarROI([1 2]) = (tmpVarROI([1 2])*params.rho_0);%+params.rho_0;
            rectangle('Position', [tmpVarROI(1) tmpVarROI(3) tmpVarROI(2)-tmpVarROI(1) tmpVarROI(4)-tmpVarROI(3)], 'EdgeColor', 'w', 'LineWidth', 1, 'LineStyle', '-.');
        end
        title(subplot_labels(tt+length(times)*2, 2));
    end
end
if plots == 1
    fig = gcf;
    figure_print_format(fig);
    set(gcf, 'PaperUnits', 'centimeters', 'Units', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 23.0452 18.0181], 'Position', [0 0 23.0452 18.0181]);
    new_dir = cd;
    if isFinal
        %figure_print_format(fig);
        cd(orig_dir);

        % IF YOUR MATLAB VERSION IS BEFORE R2020a COMMENT THESE LINES
        exportgraphics(fig, ['../02_Raster/', output_fnm, '.png'], 'Resolution', 300);
        exportgraphics(fig, ['../03_Vector/', output_fnm, '.eps']);
    end
    %pause;
    cd(new_dir)
end

if plots == 3
    figure
    usp_to_physical(times(end), 'rho', 'enstrophy', [xlims ylims], var_lims, var_ROI);
    set(gcf, 'PaperUnits', 'centimeters', 'Units', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 23.0452 18.0181], 'Position', [0 0 23.0452 18.0181]);
elseif plots ==2
    figure
    ROI = usp_to_physical(times(end), 'rho', 'enstrophy', [xlims ylims], var_lims, var_ROI);
    clf
    tlo = tiledlayout(3, 1,  'TileSpacing', 'compact', 'Padding', 'Compact');
    nexttile;
    %pcolor(x, z, (rho*params.rho_0)+params.rho_0);cmocean('dense')
    pcolor(x, z, (rho*params.rho_0));cmocean('dense')

    %var_lims([3 4]) = round((var_lims([3 4])*params.rho_0)+params.rho_0);
    var_lims([3 4]) = round((var_lims([3 4])*params.rho_0));

    var_ROI([1 2]) = var_ROI([1 2])*params.rho_0;
    colormap(gca, cmocean('dense'));
    caxis(var_lims([3 4]));
    c = colorbar('location', 'EastOutside');
    ylabel(c, "$\rho' (kg m^{-3})$"); ylabel('z (m)');
    axis([xlims ylims]);
    xticklabels([]);
    hold on
    plot(x(:, 1), z(:, 1), 'k-');

    nexttile;
    pcolor(ROI.x, ROI.z, double(~ROI.Region));
    caxis([0 1]); cmocean('amp');
    ylabel('$z (m)$'); xlabel('$x (m)$')
    hold on
    plot(x(:, 1), z(:, 1), 'k-');
    axis([xlims ylims]);

    nexttile;
    %imagesc((myvar1*params.rho_0)+params.rho_0, myvar2, log10(usp)); set(gca, 'YDir', 'normal');
    imagesc((myvar1*params.rho_0), myvar2, log10(usp)); set(gca, 'YDir', 'normal');

    colormap(gca, 'plasma'); caxis([-6 -2]);
    c = colorbar; ylabel(c, '$\log_{10}(W)$');
    axis(var_lims([3 4 1 2]));
    hold on
    rectangle('Position', [var_ROI(1) var_ROI(3) var_ROI(2)-var_ROI(1) var_ROI(4)-var_ROI(3)], 'EdgeColor', 'w', 'LineWidth', 1, 'LineStyle', '-.');
    box on; axis square;
    xticks(round([var_lims(3), mean(var_lims([3 4])), var_lims(4)]))
    ylabel('$\Omega (s^{-2})$'); xlabel("$\rho' (kg m^{-3})$")
    set(gcf, 'PaperUnits', 'centimeters', 'Units', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 11.6 11.1], 'Position', [0 0 11.6 11.1]);
end
fig = gcf;
if isFinal && ~~ plots
    fig = gcf;
    figure_print_format(fig);

    cd(orig_dir);

    % IF YOUR MATLAB VERSION IS BEFORE R2020a COMMENT THESE LINES
    exportgraphics(fig, ['../02_Raster/', output_fnm, '.png'], 'Resolution', 300);
    exportgraphics(fig, ['../03_Vector/', output_fnm, '.eps']);
    %export_fig(fig, ['../03_Vector/', output_fnm, '.eps'], '-eps');

end

