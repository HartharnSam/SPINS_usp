function data = get_spins_3Ddata(var, ii, xminInd, xmaxInd, yminInd, ymaxInd, zInds)
%GET_SPINS_3DDATA - returns 3D spins data, basically a wrapper for
%spins_reader_new, but it can do some computations for certain derived variables

switch lower(var)
    case 's'
        data = spins_reader_new('s', ii, xminInd:xmaxInd, yminInd:ymaxInd, zInds);
        data = data.*(data > 0);
    case 'enstrophy'
        try
            data = spins_reader_new('enst', ii, xminInd:xmaxInd,  yminInd:ymaxInd,  zInds);
        catch
            data = 0.5*spins_reader_new('vorty', ii, xminInd:xmaxInd, yminInd:ymaxInd,  zInds).^2;
        end
    case 'ke'
        u = spins_reader_new('u', ii, xminInd:xmaxInd,  yminInd:ymaxInd, zInds);
        w = spins_reader_new('w', ii, xminInd:xmaxInd,  yminInd:ymaxInd, zInds);
        v = spins_reader_new('v', ii, xminInd:xmaxInd,  yminInd:ymaxInd, zInds);
        data = 0.5*(u.^2 + w.^2 + v.^2);
    case 'ke_v'
        v = spins_reader_new('v', ii, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        data = 0.5.*(v.^2);
    case 'vorty'
        try
            data = spins_reader_new("vorty", ii, xminInd:xmaxInd, yminInd:ymaxInd, zInds);
        catch
            data = spins_derivs('vorty', ii);
            warning("Using spins_derivs, no save")
            data = data(xminInd:xmaxInd, yminInd:ymaxInd, zInds);
        end
    case 'speed'
        u = spins_reader_new('u',ii, xminInd:xmaxInd, yminInd:ymaxInd, zInds);
        w = spins_reader_new('w',ii, xminInd:xmaxInd, yminInd:ymaxInd, zInds);
        v = spins_reader_new('v', ii, xminInd:xmaxInd, yminInd:ymaxInd, zminInd:zmaxInd);
        data = sqrt(u.^2 + w.^2 + v.^2);
    case 'diss'
        data = spins_reader_new('diss', ii, xminInd:xmaxInd, yminInd:ymaxInd, zInds);
        data = log10(data);
    case 'rho'
        try
            data = spins_reader_new('rho', ii, xminInd:xmaxInd, yminInd:ymaxInd, zInds);
        catch
            %rho0 = params.rho_0;
            data = (eqn_of_state(spins_reader_new('t', ii, xminInd:xmaxInd,...
                yminInd:ymaxInd, zInds), 0));
        end
        %data1 = data1 + -0.5*params.delta_rho * tanh((z-params.rho_loc)/params.dz_rho);
    case 'rho_z2'
        try
            data = spins_reader_new('rho_z', ii, xminInd:xmaxInd, yminInd:ymaxInd, zInds).^2;
        catch
            spins_derivs('rho_z', ii, true);
            data = spins_reader_new('rho_z', ii, xminInd:xmaxInd, yminInd:ymaxInd, zInds).^2;
        end
    case 'rho_zz2'
        try
            data = log10(spins_reader_new('rho_zz', ii, xminInd:xmaxInd, yminInd:ymaxInd, zInds).^2);
        catch
            spins_derivs('rho_zz', ii, true);
            data = log10(spins_reader_new('rho_zz', ii, xminInd:xmaxInd, yminInd:ymaxInd, zInds).^2);
        end
        
    otherwise
        try
            data = spins_reader_new(var, ii, xminInd:xmaxInd, yminInd:ymaxInd, zInds);
        catch
            error([var, ' not configured']);
        end
end
end