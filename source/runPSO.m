
function [fitted_params, fval, exitflag, message] = runPSO(objfun, nva)

    arguments
        nva.setup = struct
        nva.sim_time = 7200
        nva.min_peak_prom = 10
        nva.min_peak_sep  = 10
        nva.peak_number   = 11 
    end
    # !!! objective function
    [fitted_params, fval, exitflag, message] = particleswarm(objfun, nvars, lb, ub, options);
end