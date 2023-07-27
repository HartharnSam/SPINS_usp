%PLOT_ALL_USP_PHYS - Produces a series of plots for comparisons between
%simulations and timesteps of USP (including Region of Interest)

clc; clearvars; close all;

spinsstartup;

plots = 1;
isFinal = true;
%% Set up each case, and plot the cases's USP timeseries
% Set up surge
data_directory{1} = "../../../../03_ShoalingStratification/02_Raw_data/Model/three_layer_case/27_111120";
times(1, :) = [NaN 62 70 75];
var_lims(1, :) = [0 2 -.0095 .0095];
var_ROI(1, :) = [-0.0085 .0085 .5 2];
xlims(1, :) = [5.5 7]; ylims(1, :) = [-.3 0];
sptitle{1} = 'Surging';
output_fnm = 'SF1_QSP_rhoE_Evolution_Surge';

i = 1;
usp_phys_plotter(data_directory{i}, times(i, ~isnan(times(i, :))), var_lims(i, :), var_ROI(i, :), xlims(i, :), ylims(i, :), output_fnm, plots, isFinal);

% Set up collapse
data_directory{2} = "../../../../03_ShoalingStratification/02_Raw_data/Model/three_layer_case/26_091120";
times(2, :) = [55 60 70 75];
var_lims(2, :) = [0 50 -.0095 .0095];
var_ROI(2, :) = [-0.0085 .0085 25 50];
xlims(2, :) = [5.5 7]; ylims(2, :) = [-.3 0];
sptitle{2} = 'Collapse';
output_fnm = 'F3_QSP_rhoE_Evolution_Collapse';

i = 2;
usp_phys_plotter(data_directory{i}, times(i, ~isnan(times(i, :))), var_lims(i, :), var_ROI(i, :), xlims(i, :), ylims(i, :), output_fnm, plots, isFinal);


% Set up plunging
data_directory{3} = "../../../../03_ShoalingStratification/02_Raw_data/Model/three_layer_case/24_071020_2";
times(3, :) = [55 60 65 70]; 
var_lims(3, :) = [0 150 -.0095 .0095];
var_ROI(3, :) = [-0.0085 .0085 75 150];
xlims(3, :) = [5.5 7]; ylims(3, :) = [-.3 0];
sptitle{3} = 'Plunging';
output_fnm = 'SF2_QSP_rhoE_Evolution_Plunge';

i = 3;
usp_phys_plotter(data_directory{i}, times(i, ~isnan(times(i, :))), var_lims(i, :), var_ROI(i, :), xlims(i, :), ylims(i, :), output_fnm, plots, isFinal);


% Set up fission
data_directory{4} = "../../../../06_Fission_Bolus/02_Raw_data/Thin_20L_33_220721/";
times(4, :) = [100 108 120 136];
var_lims(4, :) = [0 50 -.0095 .0095];
xlims(4, :) = [9 12]; ylims(4, :) = [-.3 0];
var_ROI(4, :) = [ -0.0085 0.0085 10 50];
sptitle{4} = 'Fission';

output_fnm = 'SF3_QSP_rhoE_Evolution_Fission';
i = 4;
usp_phys_plotter(data_directory{i}, times(i, ~isnan(times(i, :))), var_lims(i, :), var_ROI(i, :), xlims(i, :), ylims(i, :), output_fnm, plots, isFinal);

%% Then plot final timestep in comparison plot
clf; close all
isFinal = true;
output_fnm = 'F4_ShoalROI';
usp_phys_plotter_comb(data_directory, times, var_lims, var_ROI, xlims, ylims, output_fnm, isFinal, sptitle);


%% Bit of script to save data
% Saving the data for archival
% cd(data_directory);
% mkdir('archival');
% casename = 'CollapsingISW';
% 
% for i = 1:length(times)
%     [usp{i}, myVar1, myVar2] = usp_2d(times(i), 'rho', 'enstrophy', [xlims ylims], var_lims);
%     time{i} = times(i);
% end
% 
% CollapsingUSP = struct('Casename', casename, 'Timestep', time, 'SpatialLimits', [xlims ylims],...
%     'var_lims', var_lims, 'rho_bins', myVar1, 'enstro_bins', myVar2,  'USP', usp);
% 
% save('archival/CollapsingUSP.mat', 'CollapsingUSP');