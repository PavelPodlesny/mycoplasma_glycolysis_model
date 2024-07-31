
function fitModelToData(model_dir, data_dir, nva)

        arguments
                model_dir char  % a directory (not a full path), where model's .xml file is stored
                data_dir  char  % a directory (not a full path), where data for model's fitting is stored
                nva.repeats_number  uint8 = 1
                nva.threads         uint8 = 1
                nva.use_parallel  logical = false
                nva.sens_analysis logical = false
                nva.pooled        logical = true % pool data together (a model will be fitted across all experimental samples)
                nva.error_model      char = 'constant'
                nva.SPECIES_UNITS    char = 'micromolarity'
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

                %% write log file
                nva.model_path = model_path;
                nva.data_path = data_path;
                nva.output_path = output_dir;
                nva.date = date;
                logger(nva);
                
                %% load model's variants (different initial conditions: 1 variant per 1 group)
                variants = getModelVariants(data_path);
                
                %% load a table with experimental data
                %  exp_data.Properties.VariableNames = {'GroupID', 'Time', ... Species ...}
                exp_data = readExperimentalData(data_path);
                assert(length(variants)==length(unique(exp_data.GroupID)), ...
                        "ERROR: number of variants does not match number of groups in exp. data!");
                %% get experimental data
                [grp_data, params2est, resp_map] = getFitData(data_path, exp_data, nva.SPECIES_UNITS);
                %% create a table where optimization's results will be saved
                opt_res_vars = {'Iteration', 'Name', 'Estimate', 'StandardError', 'Bounds'};
                opt_res_table = table('Size', [0, 5], ...
                                        'VariableTypes', {'int32', 'string', 'double', 'double', 'double'}, ...
                                        'VariableNames', opt_res_vars);
                                        
                disp('Start optimization...');
                %% optimization                        
                for (iter=1:nva.repeats_number)
                        page_is_1st = true;
                        plot_path = [output_dir filesep 'fit_exp_iter_' num2str(iter) '.pdf'];
                        
                        disp(['   ' datestr(now, 'dd-mmm-HH-MM') ': iter. ' num2str(iter) ' starts...']);
                        tic;
                        optim_results = sbiofit(model, grp_data, resp_map, params2est, [], 'scattersearch', ...
                                                'Variants', variants, 'UseParallel', nva.use_parallel, ...
                                                'Pooled', nva.pooled, 'ErrorModel', nva.error_model, ...
                                                'SensitivityAnalysis', nva.sens_analysis);
                        elapsed_time = toc;
                        disp(['      iter. ' num2str(iter) ' finished: elapsed time = ' num2str(elapsed_time) ' sec.']);
                        
                        % update a table with optimization's results
                        [sim_data, est_params] = fitted(optim_results);
                        % !!! sim_data is an array
                        %% plot graphs
                        
                        for (i=1:width(variants))
                                var = variants{i};
                                grp_var_data = exp_data(exp_data.GroupID==str2num(var.Name),:);
                                plotSimData(sim_data(i), plot_path, page_is_1st,...
                                        plot_exp_data = true, grp_data = grp_var_data, resp_map = resp_map);
                                if page_is_1st
                                        page_is_1st = false;
                                end
                        end
                        %% debug
                        %disp(est_params);
                        %disp(varfun(@class, est_params, 'OutputFormat', 'cell'));
                        %% update the table 
                        iter_col = repmat(iter, height(est_params), 1);
                        est_params.Iteration = iter_col;
                        opt_res_table = vertcat(opt_res_table, est_params);
                        save([output_dir filesep 'OptimResults_iter' num2str(iter) '.mat'], '-v7.3', 'optim_results');
                end
                                
                %end
                
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

function [grp_data, params2est, resp_map] = getFitData(data_path, exp_data, species_unit)
        %!!! IMPORTANT AGREEMENT !!!: names of variables in a model and an experiment are the same for convinience,
        % so the response map has only formal sense
        %% create estimatedInfo object
        opts = detectImportOptions([data_path filesep 'Bounds.tsv'], 'FileType', 'text', ...
                                   'TextType', 'string', 'Delimiter', '\t', ...
                                   'ReadVariableNames', true, 'ReadRowNames', false);
        tbl = readtable([data_path filesep 'Bounds.tsv'], opts);
        params = (tbl.Parameter)';
        bounds = [tbl.LowerBound, tbl.UpperBound];
        params2est = estimatedInfo(params, 'Bounds', bounds);

        %% create Response Map
        exp_vars = exp_data.Properties.VariableNames(3:end); %skip GroupID and Time variables
        resp_map = strcat(string(exp_vars), " = ", string(exp_vars));

        %% create Grouped Data object
        units = repmat({species_unit}, 1, size(exp_vars, 2));
        grp_data = groupedData(exp_data, 'GroupID', 'Time'); %grouped data
        grp_data.Properties.VariableUnits = [{'', 'second'}, units];
        %%
end


function [output_dir, date] = createOutputDirectory(data_path)
%% the function creates a directory where output data will be stored
        date = datestr(now, 'dd-mmm-yy-HH-MM-SS');
        output_dir = [data_path filesep 'fit-parameters-' date];
        if ~isdir(output_dir)
                mkdir(output_dir);
        end
end

function vars = getModelVariants(data_path)
        %% create Variant object (initial cond.) for each group
        %  variants are stored as a cell array
        opts = detectImportOptions([data_path filesep 'Variant.tsv'], 'FileType', 'text', ...
                                   'TextType', 'string', 'Delimiter', '\t', ...
                                   'ReadVariableNames', true, 'ReadRowNames', false);
        tbl = readtable([data_path filesep 'Variant.tsv'], opts);
        
        ids = unique(tbl.ID);
        for (i=1:length(ids))
                var_tbl = table2cell(tbl(tbl.ID==ids(i), 2:end));
                var_cell = cell(1, size(var_tbl, 1));
                for j=1:size(var_tbl, 1)
                        var_cell{j} = var_tbl(j,:);
                end
                % !!! rewrite this part to work with not numeric VariantID
                vars{i} = sbiovariant(num2str(ids(i)), var_cell);
        end
end

function data = readExperimentalData(data_path)
        opts = detectImportOptions([data_path filesep 'ExperimentalData.tsv'], 'FileType', 'text', ...
                                   'TextType', 'string', 'Delimiter', '\t', ...
                                   'ReadVariableNames', true, 'ReadRowNames', false);
        data = readtable([data_path filesep 'ExperimentalData.tsv'], opts);
end

function logger(log_)
        fid = fopen([log_.output_path filesep 'log.txt'], 'w');
        fprintf(fid, '### Files information ###\n\tFID = %d\n', fid);
        fprintf(fid, '\tcreation date: %s\n', log_.date);
        fprintf(fid, '\tmodel dir.   : %s\n', log_.model_path);
        fprintf(fid, '\tdata  dir.   : %s\n', log_.data_path);
        fprintf(fid, '\tresults dir. : %s\n', log_.output_path);
        fprintf(fid, '### Experiment information ###\n');
        fprintf(fid, '\titer. num.   : %d\n' , log_.repeats_number);
        fprintf(fid, '\tuse parallel : %d\n' , log_.use_parallel);
        fprintf(fid, '\tpooled data  : %d\n' , log_.pooled);
        fprintf(fid, '\tsens.analysis: %d\n' , log_.sens_analysis);
        fprintf(fid, '\terror model  : %s\n' , log_.error_model);
        fclose(fid);
end
