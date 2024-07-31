function plotSimData(sim_data, plot_path, page_is_1st, nva)

        arguments
                sim_data     % SimData object
                plot_path
                page_is_1st logical
                nva.plot_exp_data logical = false
                nva.grp_data % table
                nva.resp_map % response map (match between species' names in a model and an experimental data)
        end

        % map models' variables to their units and exp. variables
        if nva.plot_exp_data
                grp_data_table = nva.grp_data;
                %% assign values to var2exp dict. (matching between species names in a model and in data)
                var2exp = dictionary;
                for i = 1:length(nva.resp_map)
                        parts = split(nva.resp_map(i), ' = ');
                        % parts{1} -- species' name in a model -- key
                        % parts{2} -- species' name in a experiment -- value
                        var2exp(parts{1}) = parts{2};
                end
        end
        
        %% assign values to var2units dict.
        var2units = dictionary;
        for i = 1:length(sim_data.DataInfo)
                var2units(sim_data.DataInfo{i}.Name) = sim_data.DataInfo{i}.Units;
        end
        
        %% plot graphs
        [time, data, names] = getdata(sim_data);
        pages_num = floor(numel(names)/4) + sign(mod(numel(names), 4)); % 4 plots per page
        for(page = 1:pages_num)
                fig = figure('Visible','off');
                for(shift = 1:4)
                        i = (page-1)*4 + shift;
                        %i = min((page-1)*4 + shift, numel(names));
                        if ( i <= numel(names) )
                                subplot(2, 2, shift);
                                plot(time, data(:,i), '-', 'lineWidth', 1);
                                
                                if (nva.plot_exp_data && isKey(var2exp, names{i}))
                                        var_name = var2exp(names{i});
                                        hold on;
                                        scatter(grp_data_table.Time, grp_data_table{:,[var_name]}, 12, "red");
                                        hold off;
                                end
                                title(names{i}, 'FontSize', 8);
                                ylabel(['Amount, ' var2units(names{i})], 'Interpreter', 'none', 'FontSize', 8);
                                xlabel('Time, sec.', 'FontSize', 8);
                        end
                end
                if page_is_1st
                        exportgraphics(fig, plot_path, 'BackgroundColor', 'none', 'ContentType', 'vector');
                        page_is_1st = false;
                else
                        exportgraphics(fig, plot_path, 'Append', true, 'BackgroundColor', ...
                                        'none', 'ContentType', 'vector');
                end 
                close(fig);
        end
        
end