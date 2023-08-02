clearvars; close all; clc;
ii = 70;
var1 = 'rho'; var2 = 'KE';

physLims = [5.5 7 -.3 0];
varLims = [0 1e-2 -0.0095 0.0095];
% Option to invert the selected region (show outside the rectangle)
isInvert = false;

%% Load in data
% Compute the qsp data
params = spins_params;
xlims = physLims([1 2]);
zlims = physLims([3 4]);
[qsp, myVar1, myKE, varLims] = qsp_mapped(ii, var1, var2, physLims, varLims);

close all;
var2Range = varLims(2)-varLims(1);
var1Range = varLims(4)-varLims(3);

var2Lims = [-1 -1 -1 -1 linspace(varLims(1), var2Range/3, 30); 500 500 500 500 linspace(varLims(2), varLims(2)-var2Range/3, 30)];
var1Lims = [-1 -1 -1 -1 linspace(varLims(3), varLims(3)+var1Range/3, 30); 1 1 1 1 linspace(varLims(4), varLims(4)-var1Range/3, 30)];

x = xgrid_reader();
xminInd = nearest_index(x(:, 1), xlims(1));
xmaxInd = nearest_index(x(:, 1), xlims(2));
x = x(xminInd:xmaxInd, :);
z = zgrid_reader(xminInd:xmaxInd, []);

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


vid = VideoWriter('C:\Users\samha\OneDrive - Newcastle University\02_PhD_Project\07_CanadaMixing\06_Communication\QSP2Phys.mp4', 'MPEG-4');
vid.FrameRate = 1;
vid.Quality = 100;
open(vid)


for ii = 1:length(var2Lims)
    
    
    RegOfInterest = ~((data1 >= var1Lims(1, ii)) & (data1 <= var1Lims(2, ii)) & (data2 >= var2Lims(1, ii))...
        & (data2 <= var2Lims(2, ii)));
    data1_temp = data1;
    data2_temp = data2;
    data1_temp(RegOfInterest) = NaN;
    data2_temp(RegOfInterest) = NaN;
    
    clf
    ax1 = subplot(3, 1, 1);
    pcolor(ax1, x, z, data1_temp); shading flat;
    colormap(gca, cmocean('dense'));
    caxis(varLims([3 4]));
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
    caxis(varLims([1 2]));
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
    axis([varLims(3) varLims(4) varLims(1) varLims(2)]);
    axis square;
    box on
    hold on
    rectangle('Position', [var1Lims(1, ii) var2Lims(1, ii) (var1Lims(2, ii)-var1Lims(1,ii)) (var2Lims(2, ii)-var2Lims(1,ii))],...
        'EdgeColor', 'w');
    xlabel(var1); ylabel(var2);
    
    drawnow; pause(.1)
    drawnow; pause(0.1); % This seems to be needed to make sure it's done it's thing before "getting frame"
    F = getframe(gcf);
    writeVideo(vid, F);
    completion(ii, 60);
end
close(vid);
