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
