function findOscillations(directory, nva)
%% nva stands for name-value arguments

%% findOscillations(dir, name-value arguments)
%% Summary. The main function for finding of oscillations.
%% Input arguments.
%        directory -- char -- a work directory (not a full path)
    arguments
        % general parameters 
        directory    char   = 'test'
        nva.info     char   = ''     % information about an experiment
        nva.threads  uint8  = 1 
        nva.sim_time uint16 = 3600   % duration of a simulation
        nva.use_parallel logical = false  % parallel computation for optimization (both for SimFunction & PSO)
        nva.repeats_number uint8 = 10     % number of repeats
        nva.optimizer_setup char = 'base' % 'base' or 'accurate' (feature is not finished yet)
        % settings for a peak calling
        % !!! add description
        nva.min_peak_prom double = 10
        nva.min_peak_sep  double = 10
        nva.SolverType       char = 'sundials'
        nva.AbsTol         double = 1e-3
        nva.RelTol         double = 1e-3
        nva.max_stall_iter_fraction double = 0.25
        nva.AbsTolScaling logical = false
        nva.terminal_log  logical = true
        %nva.peak_number   = 11 
    end

    cleanupObj = onCleanup(@cleanMeUp);
    
    try
    %% constants 
    data_path = ['/mnt/ecell_data/data2/input' filesep directory];
    time_span = [linspace(0,10,101) (11:nva.sim_time)];
    %% create a directory for output data
    [output_dir, date] = createOutputDirectory(data_path);
    if nva.terminal_log
            %log_path = [output_dir filesep 'terminal_log.txt'];
            diary /mnt/ecell_data/data2/log/terminal_log.txt;
            diary on;
            disp(['#### ' date ' ####']);
            disp(['#### ' nva.info ' ####']);
    end
    disp('Reset sbioroot...');
    reset(sbioroot);
    %% subfunctions !!!
    % K = @(x) ( 2/abs(2-x-sqrt(x*x-4*x)) );
    %% load a model
    models_list = dir([data_path filesep '*.xml']);
    model = sbmlimport([data_path filesep models_list(1).name]);
    %% get a cell array of model's species
    species_number = height(model.Species);
    species_list = cell(1, species_number);
    for i=1:species_number
        species_list{i} = model.Species(i).Name;
    end
    %% read a list of parameters/variables (with boundaries) wich will be varied 
    [params, lb, ub] = getParametersListAndBounds(data_path);
    % NOTE: dim(params) = dim(lb) = dim(ub) = [1, nvars]
    nvars = width(params);
    
    %% create log file
    nva.params = params;
    nva.lb = lb;
    nva.ub = ub;
    nva.date = date;
    nva.dir = output_dir;
    logger('header', nva);
    %% create SimFunction object
    % default solver - ode15s
    disp('Configurate solver...');
    configset = getconfigset(model, 'default');
    configset.SolverType = nva.SolverType;
    configset.SolverOptions.AbsoluteTolerance = nva.AbsTol;
    configset.SolverOptions.RelativeTolerance = nva.RelTol;
    configset.SolverOptions.AbsoluteToleranceScaling = nva.AbsTolScaling;
    %model_configs.SolverOptions.OutputTimes = time_span;
    if nva.use_parallel
        delete(gcp("nocreate"));
        disp('Create parallel pool on cluster...');
        poolobj = parpool("Processes", nva.threads); % 'Threads' don't work correctly  
        %poolobj = parpool("Threads", nva.threads);
        % !!!doesn't worl properly with SimFunction
        % addAttachedFiles(poolobj,"/home/podlesnyi/.MathWorks/SimBiology/aa9a620f-7796-4ed6-90ac-dd9c8e861ccf/");
        % /home/podlesnyi/.MathWorks/SimBiology
    else
        nva.threads = 0;
    end
    modelfun = createSimFunction(model, params, species_list, [], 'UseParallel', nva.use_parallel);
    accelerate(modelfun); % compile to C code
    %% define an objective function handler
    objfun = @(p) objectiveFunction(p, modelfun, time_span, ...
                                    nva.min_peak_sep, nva.min_peak_prom);
    %objfun = @(p) objectiveFunction(p, modelfun, time_span, ...
    %                               nva.peak_number, nva.min_peak_sep, nva.min_peak_prom);
    %objfun = @objectivefun;
    %% setup an optimizer
    pso_settings.mode = nva.optimizer_setup;
    %pso_settings.objective_limit = 0;
    pso_settings.use_parallel = nva.use_parallel;
    
    if strcmp(nva.optimizer_setup, 'accurate')
        disp('Setting pso...');
        pso_settings.swarm_size  = max(100, 10*nvars);
        pso_settings.max_stall_iterations = nva.max_stall_iter_fraction*(200*nvars); % 200*nvars = default MaxIterations
        % settings from a article 
        %pso_settings.self_weight = 1.49;
        %pso_settings.social_weight = 1.49;
        %pso_settings.inertia = [0.1, 1.1];
        %pso_settings.min_neighbors_fraction = 0.8;
    end
    opts = setupPSO(pso_settings);
    %% start optimization
    % prepare a table, where optimization's results will be saved
    optim_results_vars = [{'#'} params {'metrics', 'exitflag'}];
    optim_results = table();
    
    %% define path for graphs
    plot_path = [output_dir filesep 'fit_plots.pdf'];
    %%% !!!
    page_is_1st = true; % flag, which indicates, that a graph on the 1st page will be plotted 
    %% optimization                        
    disp('Start optimization...');
    tot_etime_start =  tic;
    %parfor (iter=1:nva.repeats_number, min(round(nva.threads/2), nva.repeats_number) )
    for (iter=1:nva.repeats_number)
    %% !!! optimization
        disp([datestr(now, 'dd-mmm-HH-MM') ': iter. ' num2str(iter) ' starts...']);
            log_ = struct();
            log_.date = date;
            log_.dir  = output_dir;
        tic;
            [fitted_params, fval, exitflag, output] = particleswarm(objfun, nvars, lb, ub, opts);
        log_.elapsed_time = toc;
        disp(['    iter. ' num2str(iter) ' finished: elapsed time = ' num2str(log_.elapsed_time) ' sec.']);
            log_.fval     = abs(fval);
            log_.iter     = iter;
            log_.exitflag = exitflag;
            log_.evalnum  = output.funccount;
            log_.message  = output.message;
            log_.fitted_params = fitted_params;
        %% temporal message part (only for debugging) !!!
        %disp('parallel: ', modelfun.UseParallel)
        %disp('accelerated: ', modelfun.isAccelerated)
        %%
        logger('iteration', log_);
        row = table();
        row(1,:) = [{iter} num2cell(fitted_params) {abs(fval), exitflag}];
        optim_results = [optim_results; row];
        %optim_results(end+1,:) = [{iter} num2cell(fitted_params) {fval, exitflag}];
%             optim_results.iter(iter) = iter;
%             optim_results.fitted_parameters(iter, 1:nvars) = fitted_params;
%             optim_results.metrics(iter) = fval;
%             optim_results.exit_flag(iter) = exitflag;
        %% save solution with new parameters
        %[t, sol] = modelfun(fitted_params, [], [], time_span);
        %solution = array2table([t{1} sol{1}]);
        sim_data = modelfun(fitted_params, [], [], time_span);
        %% add plot
        plotSimData(sim_data, plot_path, page_is_1st);
        if page_is_1st
                page_is_1st = false;
        end 
        %close(fig);
        %% save sim. results
        [t, sol, ~] = getdata(sim_data(1));
        solution = array2table([t sol]);
        solution.Properties.VariableNames = [{'time'} species_list];
        writetable(solution, [output_dir filesep 'solution-' num2str(iter) '.csv'], ...
                   'FileType', 'text', ...
                   'Delimiter', ',');
        %%
    end
    tot_elapsed_time = toc(tot_etime_start);
    disp(['Total elapsed time: ' num2str(round(tot_elapsed_time, 0)) ' sec.']);
    disp('Job is done!');
    catch ME
        delete(gcp("nocreate"));
        diary off;
        fclose('all');
        close all;
        rethrow(ME);
    end

    function cleanMeUp()
        disp('Cleaning up...');
        disp('    Shutting down a pool of workers on a cluster...');
        delete(gcp("nocreate"));
        %% save table with results
        disp('    Saving results...');
        optim_results.Properties.VariableNames = optim_results_vars; %!!!
        writetable(optim_results, [output_dir filesep 'optimization_results.csv'], ...
                   'FileType', 'text', ...
                   'Delimiter', ',');
        disp('    Closing files...');
        diary off;
        fclose('all');
        close all;
    end
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
    opts = detectImportOptions([data_path filesep 'Bounds.tsv'], 'FileType', 'text', ...
                                   'TextType', 'char', 'Delimiter', '\t', ...
                                   'ReadVariableNames', true, 'ReadRowNames', false);
    t = readtable([data_path filesep 'Bounds.tsv'], opts);
    parameters_list = (t{:, 1})';
    lower_bounds = (t{:, 2})';
    upper_bounds = (t{:, 3})';
end

%function metrics = objectiveFunction(p, modelfun, time_span, peak_number, min_peak_sep, min_peak_prom)
function metrics = objectiveFunction(p, modelfun, time_span, min_peak_sep, min_peak_prom)
    [~, sol] = modelfun(p, [], [], time_span);
    pks = islocalmax(sol{1}, 1, 'MinProminence', min_peak_prom,...
                                'MinSeparation', min_peak_sep,...
                                'SamplePoints' , time_span);
    pks_max_number = max(sum(pks,1));
    metrics = -pks_max_number;
    
end

function opts = setupPSO(s)
    switch s.mode
        case 'accurate'
            opts = optimoptions('particleswarm', ...
                                'Display'               , 'none'                   ,...
                                'SwarmSize'             , s.swarm_size             ,...
                                'MaxStallIterations'    , s.max_stall_iterations   ,...
                                'UseParallel'           , s.use_parallel           ,...
                                'UseVectorized'         , false);
                                %'ObjectiveLimit'        , s.objective_limit        ,...
                                %'InertiaRange'          , s.inertia                ,...
                                %'SelfAdjustmentWeight'  , s.self_weight            ,...
                                %'SocialAdjustmentWeight', s.social_weight          ,...
                                %'MinNeighborsFraction'  , s.min_neighbors_fraction ,...
                                %'UseParallel'           , s.use_parallel           ,...
                                %'MaxStallIterations'    , s.max_stall_iter         ,...
                                %'MaxTime'               , s.max_time
        case 'base'
            opts = optimoptions('particleswarm', ...
                                'Display'               , 'none'                   ,...
                                'UseParallel'           , s.use_parallel           ,...
                                'UseVectorized'         , false);
                                %'ObjectiveLimit'        , s.objective_limit        ,...
                                %'UseParallel'           , s.use_parallel           ,...
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
                fprintf(fid, '\tsim. time : %d sec.\n', log_.sim_time);
                fprintf(fid, '\titer. num.: %d\n', log_.repeats_number);
                fmt = ['\tvariable parameters = [', repmat('%s, ', 1, numel(log_.params)-1), '%s]\n'];
                fprintf(fid, fmt, log_.params{:});
                fmt = ['\tupper bounds = [', repmat('%g, ', 1, numel(log_.ub)-1), '%g]\n'];
                fprintf(fid, fmt, log_.ub);
                fmt = ['\tlower bounds = [', repmat('%g, ', 1, numel(log_.lb)-1), '%g]\n'];
                fprintf(fid, fmt, log_.lb);
                fprintf(fid, '###  Solver information ###\n');
                fprintf(fid, '\tsolver type : %s\n' , log_.SolverType);
                fprintf(fid, '\tabs. tol    : %e\n' , log_.AbsTol);
                fprintf(fid, '\trel. tol    : %e\n' , log_.RelTol);
                fprintf(fid, '\tabs. tol. scaling : %d\n' , log_.AbsTolScaling);
                fprintf(fid, '###  Optimization information ###\n');
                fprintf(fid, '\tsetup type  : %s\n' , log_.optimizer_setup);
                fprintf(fid, '\tuse parallel: %d\n' , log_.use_parallel);
                fprintf(fid, '\tmin. peak sep.  : %d\n', log_.min_peak_sep);
                fprintf(fid, '\tmin. peak prom. : %d\n', log_.min_peak_prom);
                fprintf(fid, '\tMaxStallIterFraction: %.2f\n' , log_.max_stall_iter_fraction);
                %fprintf(fid, '\tpeaks number = %d\n', log_.peak_number);
            fclose(fid);
        case 'iteration'
            fid = fopen([log_.dir filesep 'log.txt'], 'a');
                date = datestr(now, 'dd-mmm-yy-HH-MM-SS');
                fprintf(fid, '%s: iter. %d finished\n', date, log_.iter);
                fprintf(fid, '\texit flag = %d\n', log_.exitflag);
                fprintf(fid, '\tmetrics   = %d\n', log_.fval);
                fprintf(fid, '\tevalnum   = %d\n', log_.evalnum);
                fprintf(fid, '\ttime      = %d min.\n', round(log_.elapsed_time/60));
                fmt = ['\tfitted parameters = [', repmat('%.2f, ', 1, numel(log_.fitted_params)-1), '%.2f]\n'];
                fprintf(fid, fmt, log_.fitted_params);
                fprintf(fid, '\tmessage: %s\n\n', log_.message);
            fclose(fid);
    end
end



%function m = objectivefun(p)
%        m = sum(p);
%end
