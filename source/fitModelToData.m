
function fitModelToData(model_dir, data_dir, nva)

        arguments
                model_dir char
                data_dir char
                nva.pooled logical = true
                nva.repeats_number = 10
                nva.use_parallel logical = false
                nva.sens_analysis = true
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
                
                models_list = dir([model_path filesep '*.xml']);
                model = sbmlimport([model_path filesep models_list(1).name]);

                [output_dir, date] = createOutputDirectory(data_path);
                [grp_data, params2est, resp_map] = getFitData(data_path, model.Species);

                disp('Start optimization...');
                optim_results = sbiofit(model, grp_data, resp_map, params2est, [], 'scattersearch', ...
                                        'Pooled', nva.pooled, 'ErrorModel', nva.error_model, ...
                                        'UseParallel', nva.use_parallel);%, 'SensitivityAnalysis', nva.sens_analysis);
                [~, est_params] = fitted(optim_results);
                
                save([output_dir filesep 'OptimResults.mat'], '-v7.3', 'optim_results');
                writetable(est_params, [output_dir filesep 'EstimatedParams.tsv'], ...
                   'FileType', 'text', ...
                   'Delimiter', '\t');
                disp('Job is gone!');
                
        catch ME
                delete(gcp("nocreate"));
                fclose('all');
                rethrow(ME);
        end

        function cleanMeUp()
                disp('Cleaning up...');
                delete(gcp("nocreate"));
                disp('    Closing files...');
                fclose('all');
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