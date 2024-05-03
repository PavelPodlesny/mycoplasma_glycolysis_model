function findOscillations(directory, nva)
%% nva stands for name-value arguments

%% findOscillations(dir, name-value arguments)
%% Summary. The main function for finding of oscillations.
%% Input arguments.
% directory   -- char      -- a work directory
% observables -- char cell -- a cell array with observed species: {'G6P', 'FDP', 'ATP'}
    dummy = 0;
    arguments
        % general parameters 
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
    objective_limit = 0;
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
    output_dir = createOutputDirectory(data_path);
    %% read a list of parameters/variables (with boundaries) wich will be varied 
    [params, lb, ub] = getParametersListAndBounds(data_path);
    nvars = width(params);
    optim_results_vars = [{'#'} params {'metrics', 'exitflag'}];
    %% !!! create log file
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
    switch nva.optimizer_setup
        case 'accurate'
            % optimizer settings
            %swarm_size  = 100;
            %self_weight = 1.49;
            %socl_weight = 1.49;
            %inertia = [0.1, 1.1];
            %min_neighbors_fraction = 0.8;
            % !!! to complete
        case 'base'
        otherwise
            % !!! handle exception
    end
    %% start optimization
    optim_results = table();
    for iter=1:nva.repeats_number
    % !!! optimization
        [fitted_params, fval, exitflag, ~] = particleswarm(objfun, nvars, lb, ub);
        optim_results(end+1,:) = [{i} num2cell(fitted_params) {fval, exitflag}];
        %% save solution with new parameters
        [t, sol] = modelfun(fitted_params, [], [], time_span);
        solution = array2table([t{1} sol{1}]);
        solution.Properties.VariableNames = [{'time'} species_list];
        writetable(solution, [output_dir filesep 'solution-' num2str(iter) '.tsv'], ...
                   'Delimiter', '\t');
        %%
    end
    optim_results.Properties.VariableNames = optim_results_vars;
    writetable(solution, [output_dir filesep 'optimization_results.tsv'], ...
               'Delimiter', '\t');
    %%plot dynamics
end

function output_dir = createOutputDirectory(data_path)
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
