function [axLab] = get_axis_labels(var, isUnits)
% GET_AXIS_LABELS - Returns an axis label for a variety of SPINS variables
% (including computed ones)
arguments
    var (1, 1) string
    isUnits (1, 1) logical = false
end
switch lower(var)
    case "rho"
        axLab = "$\rho$";
        lab_units = "$(kg m^{-3})$";
    case "rho_z"
        axLab = "$\partial\rho / \partial z$";
        lab_units = "$(kg m^{-4})$";
    case "rho_zz"
        axLab = "$\partial^2\rho / \partial z^2$";
        lab_units = "$(kg m^{-5})$";
    case "speed"
        axLab = "$|u|$";
        lab_units = "$(m s^{-1})$";
    case "diss"
        axLab = "$\epsilon$";
        lab_units = [];
    case "ke"
        axLab = "KE";
        lab_units = "$(J kg^{-1})$";
    case "enstrophy"
        axLab = "$\Omega$";
        lab_units = "$(s^{-2})$";
    case "vorty"
        axLab = "$\omega";
        lab_units = "$(s^{-1})$";
    otherwise
        axLab = var;
        lab_units = "";
end

if isUnits
    axLab = axLab+lab_units;
    axLab = replace(axLab, "$$"," ");
end
