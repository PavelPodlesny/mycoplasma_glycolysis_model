
function [bioreps, bioreps_sd] = plotModelCurve(nva)

        arguments
                nva.fit_res_path      char
                nva.exp_data_path     char
                nva.plot_path         char
                nva.rct_ID            char
                nva.rct_ID_m          char
                nva.conc_unit         char = 'mM'
                nva.ci_type           char = 'gaussian'
                nva.plot_model_dat logical = true
                nva.alpha           double =  0.1
                nva.bt_smpls_num    uint16 =  500
        end
        
        %% load model's prediction, calculate CI and scale data if is's needed
        [model_data, flag] = getFitResults(nva);
        %% load experimental data
        data = getExpData(nva);
        met_names = string(data.Properties.VariableNames); %!!!
        met_names = met_names(3:end);
        rct_IDs_e = split(data.GroupID, '_');
        % select data for the desired reaction
        rct_filter = (rct_IDs_e(:,1) == nva.rct_ID);
        
        if (width(rct_IDs_e) == 1)
                % add dummy biorep
                rct_IDs_e = [rct_IDs_e, repmat("1", [height(rct_IDs_e) 1])];  
        end
        biorep_IDs = unique(rct_IDs_e(rct_filter, 2));
        bioreps = cell(1, length(biorep_IDs));
        bioreps_sd = cell(1, length(biorep_IDs));

        for i=1:length(biorep_IDs)
                biorep_filter = (rct_IDs_e(:,2) == biorep_IDs(i));
                if (width(rct_IDs_e) == 3)
                        techrep_IDs = unique(rct_IDs_e(rct_filter & biorep_filter,3));
                        techreps = cell(1, length(techrep_IDs));
                        for j=1:length(techrep_IDs)
                                techrep_filter = (rct_IDs_e(:,3) == techrep_IDs(j));
                                techreps{j} = data{rct_filter & biorep_filter & techrep_filter, 3:end};
                        end
                        techreps_conc = cat(3, techreps{:});
                        techreps_mean = mean(techreps_conc, 3);
                        time = data.Time(rct_filter & biorep_filter & techrep_filter);
                        bioreps{i} = [time techreps_mean];
                        bioreps_sd{i} = std(techreps_conc, 0, 3);
                        %%%% depricated
                        %if (length(techrep_IDs)==1)
                        %        bioreps_sd{i} = double.empty;
                        %else
                        %        bioreps_sd{i} = std(techreps_conc, 0, 3);
                        %end
                        %%%%
                else
                        bioreps{i} = data{rct_filter & biorep_filter, 2:end};
                        bioreps_sd{i} = zeros(nnz(rct_filter & biorep_filter), width(data)-2);
                end
                
        end
        %%% get init cond.
        %initc = zeros(1, length(met_names));
        %%% TEMPORAL SOLUTION, UPDATE IS NEEDED!!!
        %met_names_ = string(model_data.Properties.VariableNames); %!!!
        %met_names_ = met_names_(3:end);
        %legend_labs = cell(size(met_names_));
        %for ii = 1:length(met_names_) 
         %     initc = round(model_data{model_data.Time==0, met_names_(ii)}, 1);
         %     legend_labs{ii} = ['[' met_names_{ii} '] = ' num2str(initc) ' ' nva.conc_unit];
        %end

        %% Plotting
        fig = figure('Visible','off');
        line_spec = {"or", "squarer", "diamondr", "^"};
        nrows = ceil(numel(met_names)/2);
        for i=1:numel(met_names)
                subplot(nrows, 2, i);
                pbaspect([1.5 1 1]);
                hold on;
                for j=1:length(bioreps)
                        c = bioreps{j};
                        sd = bioreps_sd{j};
                        e = errorbar(c(:,1), c(:,i+1), sd(:,i), line_spec{j}, 'MarkerSize', 6);
                        e.CapSize = 0;
                        axis padded;
                end
                
                if nva.plot_model_dat
                        switch flag
                                case "mat"
                                        model_data_mets = split(string(model_data.Response),'.');
                                        c = model_data(model_data_mets(:,2) == met_names(i), :);
                                        %disp([c.Time c.ConfidenceInterval(:,2) c.ConfidenceInterval(:,1)])
                                        fill([c.Time; flipud(c.Time)], ...
                                             [c.ConfidenceInterval(:,1); flipud(c.ConfidenceInterval(:,2))], ...
                                             'blue', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                                        plot(c.Time, c.Estimate, '-b', 'LineWidth', 2);
                                case "csv"
                                        plot(model_data.Time, model_data{:, met_names(i)}, '-b', 'LineWidth', 2);
                        end
                end
                dummy = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
                %legend_text = strjoin(legend_labs, '\n');
                %lgd = legend([dummy], legend_text, 'Location', 'best', 'FontSize', 7);
                %title(lgd,'Initial conc.:');
                %legend('boxoff');
                %%%
                ax = gca;
                ax.FontSize = 8;
                xlabel('Time, min.', 'FontSize', 10);
                ylabel(strcat("[", met_names(i), "], ", nva.conc_unit), 'Interpreter', 'none', 'FontSize', 10);
                hold off;
        end
        exportgraphics(fig, nva.plot_path, 'BackgroundColor', 'none', 'ContentType', 'vector');
end

function [model_data, flag] = getFitResults(nva)

        scale_factor = dictionary(["mM", "mkM"], ...
                                   [1e3, 1]);
        
        tmp = split(nva.fit_res_path, '/');
        tmp = split(tmp(end), '.');
        flag = tmp{2};

        switch flag
                case "mat"
                        %% place holder
                        fit_res = load(nva.fit_res_path).optim_results;
                        ci = sbiopredictionci(fit_res, 'Type', nva.ci_type, 'Alpha', nva.alpha, ...
                                'NumSamples', nva.bt_smpls_num);
                        model_data = ci.Results(ci.Results.Group == nva.rct_ID_m, :);
                        % scale conc.:
                        model_data{:, 4:end} = model_data{:, 4:end} ./ scale_factor(nva.conc_unit);
                case "csv"
                        opts = detectImportOptions(nva.fit_res_path, ...
                                   'FileType', 'text', ...
                                   'TextType', 'string', ...
                                   'Delimiter', ',', ...
                                   'ReadVariableNames', true, 'ReadRowNames', false);
                        opts.VariableTypes = [{'string'}, repmat({'double'}, 1, length(opts.VariableNames) - 1)];
                        fit_res = readtable(nva.fit_res_path, opts);
                        model_data = fit_res(fit_res.GroupID == nva.rct_ID_m, :);
                        % scale conc.:
                        model_data{:, 3:end} = model_data{:, 3:end} ./ scale_factor(nva.conc_unit);
                        
                otherwise
                        error('plotModelCurve:incorrect file Type', ...
                                "Error. \nFile with model's estimate should be .csv or .mat, not .%s", flag);
        end
        % scale time: sec --> min
        model_data.Time = model_data.Time ./ 60;
end

function data = getExpData(nva)
        scale_factor = dictionary(["mM", "mkM"], ...
                                   [1e3, 1]);
        
        opts = detectImportOptions(nva.exp_data_path, ...
                                   'FileType', 'text', ...
                                   'TextType', 'string', ...
                                   'Delimiter', '\t', ...
                                   'ReadVariableNames', true, 'ReadRowNames', false);
        opts.VariableTypes = [{'string'}, repmat({'double'}, 1, length(opts.VariableNames) - 1)];
        data = readtable(nva.exp_data_path, opts);

        data.Time = data.Time ./ 60; % sec. --> min.
        data{:, 3:end} = data{:, 3:end} ./ scale_factor(nva.conc_unit);
end