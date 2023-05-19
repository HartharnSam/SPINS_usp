function make_qsp_isos(var1, var2, ii, var1lim, var2lim, savefnm)
%MAKE_QSP_ISOS - Makes isosurfaces from qsp physical region.
% Nice qualitative qsp_to_physical but for 3D, which can then be used with
% plot_qsp_iso. Ideall this part is run on the remote cluster, and
% plot_qsp_iso is run locally
%
% Inputs:
%    input1 - Description
%    input2 - Description
%    input3 - Description
%
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Other m-files required: spins_reader_new
% Subfunctions: none
% MAT-files required: none
%
% See also: PLOT_QSP_ISO
% Author: Sam Hartharn-Evans
% School of Mathematics, Statistics and Physics, Newcastle University
% email address: s.hartharn-evans2@newcastle.ac.uk
% GitHub: https://github.com/HartharnSam
% 19-May-2023; Last revision: 19-May-2023
% MATLAB Version: 9.12.0.2170939 (R2022a) Update 6

%---------------------------------------------------
%% BEGIN CODE %%
%---------------------------------------------------
%clc; clearvars; close all;

var1data = spins_reader_new(var1, ii);
switch lower(var2)
    case 'ke_v'
        var2data = spins_reader_new('v', ii);
        var2data = 0.5*var2data.*var2data;
    otherwise
        var2data = spins_reader_new(var2, ii);
end

RegOfInterest = (var1data>var1lim(1) & var1data<var1lim(2)) & var2data > var2lim(1) & var2data < var2lim(2);

save(savefnm, 'RegOfInterest');


%---------------------------------------------------
%% END OF CODE %%
% --------------------------------------------------