function findOscillations(directory, nva)
%% nva stands for name-value arguments

%% findOscillations(dir, name-value arguments)
%% Summary. The main function for finding of oscillations.
%% Input arguments.
% directory   -- char      -- a work directory
% observables -- char cell -- a cell array with observed species: {'G6P', 'FDP', 'ATP'}
    arguments
        % general parameters 
        directory char
        nva.info            = ''     % information about an experiment 
        nva.sim_time        = 7200   % duration of a simulation
        nva.use_parallel    = false  % parallel computation for optimization
        nva.repeats_number  = 1      % number of repeats
        nva.optimizer_setup = 'base' % 'base' or 'accurate'
        % settings for a peak calling
        % !!! add description
        nva.min_peak_prom = 10
        nva.min_peak_sep  = 10
        nva.peak_number   = 11 
    end

    
    %% subfunctions
    % K = @(x) ( 2/abs(2-x-sqrt(x*x-4*x)) );
    %% constants 
    data_path = ['/mnt/ecell_data/data2/input' filesep directory];
    time_span = [linspace(0,10,101) (11:nva.sim_time)];
    %% load a model
    models_list = dir([data_path filesep '*.xml']);
    model = sbmlimport([data_path filesep models_list(1).name]);
    %%get a cell array of model's species
    species_number = height(model.Species);
    species_list = cell(1, species_number);
    for i=1:species_number
        species_list{i} = model.Species(i).Name;
    end
    %% create a directory for output data
    [output_dir, date] = createOutputDirectory(data_path);
    log_.date = date;
    log_.dir = output_dir;
    %% read a list of parameters/variables (with boundaries) wich will be varied 
    [params, lb, ub] = getParametersListAndBounds(data_path);
    nvars = width(params);
    
    %% !!! create log file
    nva.params = params;
    nva.lb = lb;
    nva.ub = ub;
    nva.date = date;
    nva.dir = output_dir;
    logger('header', nva);
    %% create SimFunction object
    % default solver - ode15s
    %model_configs = getconfigset(model);
    %model_configs.StopTime = nva.sim_time; #set duration of a simulation
    %model_configs.SolverOptions.OutputTimes = time_span;
    modelfun = createSimFunction(model, params, species_list, []);
    %% define an objective function handler
    objfun = @(p) objectiveFunction(p, modelfun, time_span, ...
                                    nva.peak_number, nva.min_peak_sep, nva.min_peak_prom);
    %% setup an optimizer
    pso_settings.mode = nva.optimizer_setup;
    pso_settings.objective_limit = 0;
    pso_settings.use_parallel = nva.use_parallel;
    
    if strcmp(nva.optimizer_setup, 'accurate')
        pso_settings.swarm_size  = 100;
        pso_settings.self_weight = 1.49;
        pso_settings.social_weight = 1.49;
        pso_settings.inertia = [0.1, 1.1];
        pso_settings.min_neighbors_fraction = 0.8;
    end
    opts = setupPSO(pso_settings);
    %% start optimization
    optim_results_vars = [{'#'} params {'metrics', 'exitflag'}];
    optim_results = table();
    
    for iter=1:nva.repeats_number
    %% !!! optimization
        tic;
            [fitted_params, fval, exitflag, output] = particleswarm(objfun, nvars, lb, ub, opts);
        log_.elapsed_time = toc;
        
            log_.iter = iter;
            log_.exitflag = exitflag;
            log_.fval = fval;
            log_.evalnum = output.funccount;
            
        logger('iteration', log_);
        optim_results(end+1,:) = [{iter} num2cell(fitted_params) {fval, exitflag}];
%             optim_results.iter(iter) = iter;
%             optim_results.fitted_parameters(iter, 1:nvars) = fitted_params;
%             optim_results.metrics(iter) = fval;
%             optim_results.exit_flag(iter) = exitflag;
        %% save solution with new parameters
        [t, sol] = modelfun(fitted_params, [], [], time_span);
        solution = array2table([t{1} sol{1}]);
        solution.Properties.VariableNames = [{'time'} species_list];
        writetable(solution, [output_dir filesep 'solution-' num2str(iter) '.tsv'], ...
                   'FileType', 'text', ...
                   'Delimiter', '\t');
    end
    %% save table with results
    optim_results.Properties.VariableNames = optim_results_vars;
    writetable(optim_results, [output_dir filesep 'optimization_results.tsv'], ...
               'FileType', 'text', ...
               'Delimiter', '\t');
    %% !!! plot dynamics
end

%% support functions

function [output_dir, date] = createOutputDirectory(data_path)
%% the function creates a directory where output data will be stored
    date = datestr(now, 'dd-mmm-yy-HH-MM-SS');
    output_dir = [data_path filesep 'oscill-search-' date];
    if ~isdir(output_dir)
        mkdir(output_dir);
    end
end

function [parameters_list, lower_bounds, upper_bounds] = getParametersListAndBounds(data_path)
    opts = detectImportOptions([data_path filesep 'bounds.tsv'], 'FileType', 'text', ...
                                   'TextType', 'char', 'Delimiter', '\t', ...
                                   'ReadVariableNames', true, 'ReadRowNames', false);
    t = readtable([data_path filesep 'bounds.tsv'], opts);
    parameters_list = (t{:, 1})';
    lower_bounds = (t{:, 2})';
    upper_bounds = (t{:, 3})';
end

function metrics = objectiveFunction(p, modelfun, time_span, peak_number, min_peak_sep, min_peak_prom)
    [~, sol] = modelfun(p, [], [], time_span);
    pks = islocalmax(sol{1}, 1, 'MinProminence', min_peak_prom,...
                                'MinSeparation', min_peak_sep,...
                                'SamplePoints' , time_span);
    pks_max_number = max(sum(pks,1));
    metrics = (pks_max_number - peak_number)^2;
    
end

function opts = setupPSO(s)
    switch s.mode
        case 'accurate'
            opts = optimoptions('particleswarm', ...
                                'SwarmSize'             , s.swarm_size             ,...
                                'InertiaRange'          , s.inertia                ,...
                                'SelfAdjustmentWeight'  , s.self_weight            ,...
                                'SocialAdjustmentWeight', s.social_weight          ,...
                                'MinNeighborsFraction'  , s.min_neighbors_fraction ,...
                                'UseParallel'           , s.use_parallel           ,...
                                'ObjectiveLimit'        , s.objective_limit);
                                %'MaxStallIterations'   , s.max_stall_iter         ,...
                                %'MaxTime'              , s.max_time
        case 'base'
            opts = optimoptions('particleswarm', ...
                                'UseParallel'           , s.use_parallel           ,...
                                'ObjectiveLimit'        , s.objective_limit);
    end
end

function logger(mode, log_)
    switch mode
        case 'header'
            fid = fopen([log_.dir filesep 'log.txt'], 'w');
                fprintf(fid, '### File information ###\n\tFID = %d\n', fid);
                fprintf(fid, '\tcreation date: %s\n', log_.date);
                fprintf(fid, '\tdirectory    : %s\n', log_.dir);
                fprintf(fid, '### Experiment information ###\n');
                fprintf(fid, ['\t%s\n'], log_.info );
                fprintf(fid, '\titer. num.: %d\n', log_.repeats_number);
                fmt = ['\tvariable parameters = [', repmat('%s, ', 1, numel(log_.params)-1), '%s]\n'];
                fprintf(fid, fmt, log_.params{:});
                fmt = ['\tupper bounds = [', repmat('%g, ', 1, numel(log_.ub)-1), '%g]\n'];
                fprintf(fid, fmt, log_.ub);
                fmt = ['\tlower bounds = [', repmat('%g, ', 1, numel(log_.lb)-1), '%g]\n'];
                fprintf(fid, fmt, log_.lb);
                fprintf(fid, '\tsetup type  : %s\n' , log_.optimizer_setup);
                fprintf(fid, '\tuse parallel: %d\n' , log_.use_parallel);
                fprintf(fid, '\tsim. time    = %d\n', log_.sim_time);
                fprintf(fid, '\tpeaks number = %d\n', log_.peak_number);
                fprintf(fid, '\tminPeakSep   = %d\n', log_.min_peak_sep);
                fprintf(fid, '\tminPeakProm  = %d\n', log_.min_peak_prom);
            fclose(fid);
        case 'iteration'
            fid = fopen([log_.dir filesep 'log.txt'], 'a');
                date = datestr(now, 'dd-mmm-yy-HH-MM-SS');
                fprintf(fid, '%s: iter. %d finished\n', date, log_.iter);
                fprintf(fid, '\texit flag = %d\n', log_.exitflag);
                fprintf(fid, '\tmetrics   = %d\n', sqrt(log_.fval));
                fprintf(fid, '\tevalnum   = %d\n', log_.evalnum);
                fprintf(fid, '\ttime      = %d min.\n', round(log_.elapsed_time/60));
            fclose(fid);
    end
end
