function fig = plotSimData(sim_data, nva)

        arguments
                sim_data     % SimData object
                nva.plot_exp_data logical = false
                nva.grp_data % groupedData object
                nva.resp_map % response map (match between species' names in a model and an experimental data)
        end

        % map models' variables to their units and exp. variables
        if nva.plot_exp_data
                grp_data_table = groupedData2table(nva.grp_data);
                %% assign values to var2exp dict.
                var2exp   = dictionary;
                for i = 1:length(nva.resp_map)
                        parts = split(nva.resp_map(i), ' = ');
                        % parts{1} -- species' name in a model -- key
                        % parts{2} -- species' name in a experiment -- value
                        var2exp(parts{1}) = parts{2};
                end
                %%
        end
        %% assign values to var2units dict.
        var2units = dictionary;
        for i = 1:length(sim_data(1).DataInfo)
                var2units(sim_data(1).DataInfo{i}.Name) = sim_data(1).DataInfo{i}.Units;
        end
        %% plot graphs
        fig = figure('Visible','off');
        [time, data, names] = getdata(sim_data(1));

        sqrtnames = sqrt(numel(names));
        nrows = round(sqrtnames);
        ncols = ceil(sqrtnames);

        for(i = 1:numel(names))
                subplot(nrows, ncols, i);
                plot(time, data(:,i), '-', 'lineWidth', 1);
                
                if (nva.plot_exp_data && isKey(var2exp, names(i)))
                        var_name = var2exp(names(i));
                        exp_data = grp_data_table(:, ['GroupID', 'Time', var_name]);
                        hold on;
                        s = scatter(exp_data, 'Time', var_name, 'ColorVariable', 'GroupID');
                        s.SizeData = 10;
                        hold off;
                end
                
                title(names(i));
                ylabel(['Amount, ' var2units(names(i))], 'Interpreter', 'none', 'FontSize', 8);
                xlabel('Time, sec.', 'FontSize', 8);
        end
end