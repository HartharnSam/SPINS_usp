function [eof, pc, error] = usp_eof(t1, t2, var1, var2, t_select, numMode, spatLims, varLims, calcUSP)
%USP_EOF - Compute EOF error maps for feature identification in USP plots,
%as a tool for identifying key timesteps
%
% Syntax:  usp_eof(t1, t2, var1, var2, spatLims, varLims, t_select, numMode, loadUSP)
% %t1 = 50;t2 = 150;
% Inputs:
%    t1 - First timestep 
%    t2 - End timestep
%    var1 - variable to compare to var2 (on x axis)
%    var2 - Variable to compare to var1 (on y axis)
%    spatLims - [optional] Spatial region of physical space [xmin xmax ymin ymax zmin zmax]
%       optionally, only set the x limits. Defaults to full size of tank
%    varLims - [optional] realistic limits of the variables to investigate as [var2min var2max var1min var2max]
%    region - [optional] USP region of interest to display data for
%
% Other m-files required: usp_2d, spins_params, cmocean, plasma, eof_error_map,
% subaxis
%
% See also: usp_2d, eof_error_map,
% Author: Sam Hartharn-Evans
% Based on the method of Shaw & Stastna (2019) -
% https://doi.org/10.1371/journal.pone.0225439, 
% code available: 
% https://git.uwaterloo.ca/j9shaw/PLOS-one-2019
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 17-Jun-2022; Last revision: 17-Aug-2023
% MATLAB Version: 9.12.0.1956245 (R2022a) Update 2
%
%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------

%% Set up optional input parameters
params = spins_params;

if nargin < 9
    if isfile('all_qsp.mat')
        calcUSP = false;
    else
        calcUSP = true;
    end
end
if nargin < 7 || isempty(spatLims)
    spatLims = [params.min_x params.Lx params.min_z params.Ly];
end

%% Load in data
if nargin < 8 || (isempty(varLims))
    [usp, myVar1, myVar2, varLims] = usp_2d(mean(t1, t2), var1, var2, spatLims);%[0 params.Lx params.min_z params.min_z + params.Lz]);
else
    [usp, myVar1, myVar2, ~] = usp_2d(mean(t1, t2), var1, var2, spatLims, varLims);%[0 params.Lx params.min_z params.min_z + params.Lz]);
end

if calcUSP
    usp_timeseries = repmat(usp.*0, [1 1 length(t1:t2)]);
    
    for ii = t1:t2
        usp_timeseries(:, :, ii-t1+1) = usp_2d(ii, var1, var2, spatLims, varLims);
        completion(ii-t1, t2-t1); % optional script to check completion status, comment if not available
    end
    save('all_qsp.mat', 'usp_timeseries');
end
load('all_qsp', 'usp_timeseries');

%% Use the SPINS eof_error_map
% This is a plot of "how many modes are needed at a given timestep to 
% re-construct the USP", or some measure of complexity

[m, n, o] = size(usp_timeseries);
error_contour = 0.05;
eof_input_data = reshape(usp_timeseries, m*n, o);
[eof_error, u, coeff] = eof_error_map(eof_input_data, o);
tiledlayout(2,1);
nexttile;
imagesc(t1:t2, 1:size(usp_timeseries,3), eof_error); shading flat; colorbar;
hold on
[~, cc] = contour(t1:t2, 1:size(usp_timeseries, 3), eof_error, [1 1].*error_contour, 'w');
colormap(plasma);
axis tight
ylim([1 20]) % Plot first 20 modes only
xlabel('$time (s)$');
set(gca, 'YDir', 'normal');
ylabel('Mode')
legend(cc, ['Error = ', num2str(error_contour)], 'AutoUpdate','off', 'Color', 'k', 'TextColor', 'w');
xline(t_select, '--w')

% Plot coefficient timeseries for feature identification
nexttile;
p1 = plot(t1:t2, coeff(1:numMode, :));
for ii = 1:numMode
    leg_labels{ii} = ['Mode-', num2str(ii)];
end
legend(p1, leg_labels)
ylabel('Coefficient');
xlabel('$time (s)$');

figure_print_format(gcf);

%% Plot the modes themselves ( generally not useful)
figure;
for ii = 1:numMode
    u2 = reshape(u(:, ii), n, m);
    subplot(numMode, 1, ii);
    temp_data = u2+(1.5*min(u2));
    temp_data(temp_data<=0) = NaN;
    imagesc(myVar1, myVar2, log10(temp_data));
    set(gca, 'YDir', 'normal');
    shading flat; axis tight;
    axis square
    colorbar;
    ylabel('KE'); xlabel('rho')
    cmocean('curl');
end

%% Plot a selected time step (t_select)
figure
subaxis(1, 2, 1)
pcolor(log10(usp_timeseries(:, :, t_select))); shading flat; colormap(plasma); colorbar;
axis tight;%caxis([0 0.05])

%% Reconstruct the data with a set number of modes
recon_temp = zeros(m*n, o); % loop variable for reconstruction

for ii = 1:numMode
    recon_temp = recon_temp + (u(:, ii)*coeff(ii, :));
end

recon_data = bsxfun(@plus, recon_temp(:, t_select), mean(eof_input_data,2)); % remove mean
recon_data(recon_data<=0) = NaN;
recon_data = reshape(recon_data, m, n);

% And plot this
subaxis(1, 2, 2)
pcolor(log10(recon_data)); shading flat; colormap(plasma)
colorbar;axis tight; caxis([-6 -.9])