
function fitModelToData(model_dir, data_dir, nva)

        arguments
                model_dir char % a directory (not a full path), where model's .xml file is stored
                data_dir char  % a directory (not a full path), where data for model's fitting is stored
                nva.pooled logical = true % pool data together (a model will be fitted across all experimental samples)
                nva.repeats_number = 10
                nva.use_parallel logical = false
                nva.sens_analysis = false
                nva.threads double = 1
                nva.error_model char = 'constant'
        end

        cleanupObj = onCleanup(@cleanMeUp);
        model_path = ['/mnt/ecell_data/data2/input' filesep model_dir];
        data_path = [model_path filesep data_dir];
        
        try
                disp('Reset sbioroot...');
                reset(sbioroot);

                if nva.use_parallel
                        disp('Create parallel pool on cluster...');
                        delete(gcp("nocreate"));
                        poolobj = parpool("Processes", nva.threads); %!!!
                else
                        nva.threads = 0;
                end

                %% load and accelerate a model 
                models_list = dir([model_path filesep '*.xml']);
                model = sbmlimport([model_path filesep models_list(1).name]);
                sbioaccelerate(model, getconfigset(model, 'default'));
                %% create output directory and read experimental data
                [output_dir, date] = createOutputDirectory(data_path);
                [grp_data, params2est, resp_map] = getFitData(data_path, model.Species);

                %nva.repeats_number*length(params2est)
                % calculate size of a output table with estimated parameters
                %% create a table where optimization's results will be saved
                opt_res_vars = {'Iteration', 'Name', 'Estimate', 'StandardError', 'Bounds'};
                opt_res_table = table('Size', [0, 5], ...
                                        'VariableTypes', {'int32', 'string', 'double', 'double', 'double'}, ...
                                        'VariableNames', opt_res_vars);
                %% define path for graphs
                plot_path = [output_dir filesep 'fit_exp_plots.pdf'];
                page_is_1st = true; % flag, which indicates, that a graph on the 1st page will be plotted 
                %% optimization                        
                disp('Start optimization...');
                for (iter=1:nva.repeats_number)
                        disp([datestr(now, 'dd-mmm-HH-MM') ': iter. ' num2str(iter) ' starts...']);
                        tic;
                        optim_results = sbiofit(model, grp_data, resp_map, params2est, [], 'scattersearch', ...
                                                'Pooled', nva.pooled, 'ErrorModel', nva.error_model, ...
                                                'UseParallel', nva.use_parallel, 'SensitivityAnalysis', nva.sens_analysis);
                        elapsed_time = toc;
                        disp(['    iter. ' num2str(iter) ' finished: elapsed time = ' num2str(elapsed_time) ' sec.']);
                        % update a table with optimization's results
                        [sim_data, est_params] = fitted(optim_results);
                        
                        %% plot graphs
                        fig = plotSimData(sim_data, plot_exp_data = true, grp_data = grp_data, resp_map = resp_map);
                        %plot(optim_results);
                        if page_is_1st
                                exportgraphics(fig, plot_path, 'BackgroundColor', 'none', 'ContentType', 'vector');
                                page_is_1st = false;
                        else
                                exportgraphics(fig, plot_path, 'Append', true, ...
                                                'BackgroundColor', 'none', 'ContentType', 'vector');
                        end 
                        close(fig);
                        %% debug
                        %disp(est_params);
                        %disp(varfun(@class, est_params, 'OutputFormat', 'cell'));
                        %% update the table 
                        iter_col = repmat(iter, height(est_params), 1);
                        est_params.Iteration = iter_col;
                        opt_res_table = vertcat(opt_res_table, est_params);
                        save([output_dir filesep 'OptimResults_iter' num2str(iter) '.mat'], '-v7.3', 'optim_results');
                end
                %writetable(est_params, [output_dir filesep 'EstimatedParams.tsv'], ...
                %   'FileType', 'text', ...
                %   'Delimiter', '\t');
                disp('Job is done!');
                
        catch ME
                delete(gcp("nocreate"));
                fclose('all');
                close all;
                rethrow(ME);
        end

        function cleanMeUp()
                disp('Cleaning up...');
                disp('    Shutting down a pool of workers on a cluster...');
                delete(gcp("nocreate"));
                disp('    Saving results...');
                writetable(opt_res_table, [output_dir filesep 'optimization_results.csv'], ...
                        'FileType', 'text', ...
                        'Delimiter', ',');
                disp('    Closing files...');
                fclose('all');
                close all;
        end
end

%% support function

function [grp_data, params2est, resp_map] = getFitData(data_path, model_species)
        %% create estimatedInfo object
        opts = detectImportOptions([data_path filesep 'Bounds.tsv'], 'FileType', 'text', ...
                                   'TextType', 'string', 'Delimiter', '\t', ...
                                   'ReadVariableNames', true, 'ReadRowNames', false);
        tbl = readtable([data_path filesep 'Bounds.tsv'], opts);
        params = (tbl.parameter)';
        bounds = [tbl.lower_bound, tbl.upper_bound];
        params2est = estimatedInfo(params, 'Bounds', bounds);
        %% create ResponseMap
        % exp_vars = (tbl.experiment_variable)';
        % resp_map = strcat(mod_vars, " = ", exp_vars);
        % 
        %% create groupedData object and ResponseMap
        %units = cellstr(["", "second", (tbl.unit)']); % ID, Time, Variables
        opts = detectImportOptions([data_path filesep 'ExperimentalData.tsv'], 'FileType', 'text', ...
                                   'TextType', 'string', 'Delimiter', '\t', ...
                                   'ReadVariableNames', true, 'ReadRowNames', false);
        tbl = readtable([data_path filesep 'ExperimentalData.tsv'], opts);
        exp_vars = tbl.Properties.VariableNames(3:end);
        
        n = height(model_species);
        mod_vars = cell(1, n);
        units = cell(1, n);
        for i=1:n
                mod_vars{i} = model_species(i).Name;
                units{i} = model_species(i).Units;
        end
        
        resp_map = strcat(string(mod_vars), " = ", string(exp_vars));
        grp_data = groupedData(tbl, 'ID', 'Time'); %grouped data
        grp_data.Properties.VariableUnits = [{'', 'second'}, units];
end


function [output_dir, date] = createOutputDirectory(data_path)
%% the function creates a directory where output data will be stored
        date = datestr(now, 'dd-mmm-yy-HH-MM-SS');
        output_dir = [data_path filesep 'fit-parameters-' date];
        if ~isdir(output_dir)
                mkdir(output_dir);
        end
end
