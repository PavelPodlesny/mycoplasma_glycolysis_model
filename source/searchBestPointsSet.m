
function searchBestPointsSet(model_dir, data_dir, max_points2skip, nva)
        arguments
        % fitModelToData arguments
                model_dir char
                data_dir  char
                max_points2skip     uint8 = 1 % searchBestPointsSet argument
                nva.error_model      char = 'constant'
                nva.SPECIES_UNITS    char = 'micromolarity'
                nva.ci_type          char = 'bootstrap'
                nva.optim_func       char = 'scattersearch'
                nva.repeats_number  uint8 = 1
                nva.threads         uint8 = 1
                nva.bt_samples_num uint16 = 500 % number of samples for bootstrapping
                nva.use_parallel_fit  logical = false
                nva.use_parallel_ci   logical = false
                nva.sens_analysis     logical = false
                nva.pooled            logical = true
        end

        %% read file with all model's variants
        model_path = ['/mnt/ecell_data/data2/input' filesep model_dir];
        data_path  = [model_path filesep data_dir];
        opts = detectImportOptions([data_path filesep 'Variant_all.tsv'], 'FileType', 'text', ...
                                   'TextType', 'string', 'Delimiter', '\t', ...
                                   'ReadVariableNames', true, 'ReadRowNames', false);
        opts.VariableTypes = {'string', 'string', 'string', 'string', 'double'};
        variants = readtable([data_path filesep 'Variant_all.tsv'], opts);
        %% read file with all model's exp. data
        opts = detectImportOptions([data_path filesep 'ExperimentalData_all.tsv'], 'FileType', 'text', ...
                                   'TextType', 'string', 'Delimiter', '\t', ...
                                   'ReadVariableNames', true, 'ReadRowNames', false);
        opts.VariableTypes = [{'string'}, repmat({'double'}, 1, length(opts.VariableNames) - 1)];
        exp_data_all = readtable([data_path filesep 'ExperimentalData_all.tsv'], opts);
        %%
        id = unique(variants.ID);
        n_id = length(id);
        %% prepare table for results
        res_vars = {'NumPoints', 'IDs', 'Dir', 'LogLikelihood', 'AIC', 'BIC', 'SSE'};
        res_table = table('Size', [0, 7], ...
                        'VariableTypes', {'uint32', 'string', 'string', 'double', 'double', 'double', 'double'}, ...
                        'VariableNames', res_vars);
        %%
        for (k=n_id:-1:(n_id-max_points2skip))
                disp(['searchBestPointsSet: choose ' num2str(k) ' points from ' num2str(n_id)]);
                sub_id = nchoosek(id, k);
                for(j=1:height(sub_id))
                        disp(['VariantID set: ' char(join(sub_id(j, :), ' '))]);
                        variant = variants(ismember(variants.ID, sub_id(j, :)), :);
                        exp_data = exp_data_all(ismember(exp_data_all.GroupID, sub_id(j, :)), :);
                        writetable(exp_data, [data_path filesep 'ExperimentalData.tsv'], ...
                                'FileType', 'text', 'Delimiter', '\t', ...
                                'WriteMode','overwrite', 'QuoteStrings', 'none');
                        writetable(variant, [data_path filesep 'Variant.tsv'], ...
                                'FileType', 'text', 'Delimiter', '\t', ...
                                'WriteMode','overwrite', 'QuoteStrings', 'none');
                        [R, metrics, path] = fitModelToData(model_dir, data_dir, ...
                                        optim_func=nva.optim_func, error_model=nva.error_model, ci_type=nva.ci_type, ...
                                        repeats_number=nva.repeats_number, threads=nva.threads, bt_samples_num=nva.bt_samples_num, ...
                                        use_parallel_fit=nva.use_parallel_fit, sens_analysis=nva.sens_analysis, pooled=nva.pooled, ...
                                        SPECIES_UNITS=nva.SPECIES_UNITS);
                        splitted_path = split(path, '/');
                        dir = string(splitted_path{end});
                        res_row = table(k, join(sub_id(j, :), '|'), dir, ...
                                        metrics("LogLikelihood"), metrics("AIC"), ...
                                        metrics("BIC"), metrics("SSE"), 'VariableNames', res_vars);
                        res_table = vertcat(res_table, res_row);
                end
        end
        writetable(res_table, [data_path filesep 'variant_optim_set.csv'], 'FileType', 'text', 'Delimiter', ',');
end