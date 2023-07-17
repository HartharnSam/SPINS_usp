

var1 = 'rho', var2 = 'rho';
params = spins_params;
spat_lims = [params.min_x params.min_x+params.Lx params.min_z params.min_z+params.Lz];
var_lims = [-0.0095 .0095 -.0095 .0095];

t1 = 0; t2 = 200;
for ii = t1:t2
    [usp, myVar1, myVar2] = usp_2d(ii, var1, var2, spat_lims, var_lims);
    var1_ts(:, ii-t1+1) = sum(usp);
end
%%
close all; 
clc;
pcolor(t1:t2, myVar1(2:end-1)', log(var1_ts(2:end-1, :)));
colormap(inferno);
colorbar;
