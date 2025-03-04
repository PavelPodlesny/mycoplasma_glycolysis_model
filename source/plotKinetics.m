
function plotKinetics(file_path, plot_path)
    % Reads a .csv file and plots V0 vs VariableValue for each Variable
    % Inputs:
    %   - filename: string, the name of the .csv file to read
    % Example:
    %   plot_csv_data('data.csv')

    % Read the table from the CSV file
    data = getKineticsData(file_path);
    plot_path = [plot_path filesep 'test_plot.pdf'];
    pageIs1st = true;
    % Extract unique Variables
    variables = unique(data.Variable);
    var_n = length(variables);
    disp(variables);
    disp(var_n);
    ncols = 2;
    fs = 10;
    aes = {'or', 'ob', 'om', 'og'};
    % Create a new figure
    % Loop through each unique Variable
    for i = 1:length(variables)
        var = variables{i}; % Current Variable
        % Filter rows corresponding to the current Variable
        subset = data(strcmp(data.Variable, var), :);
        fig = figure('Visible','off');
        %subplot(var_n, fix((i-1)/ncols)+1, mod((i-1),ncols)+1);
        %disp(['row: ' num2str(fix((i-1)/ncols)+1) ' col: ' num2str(mod((i-1),ncols)+1)]);
        ax = gca; ax.FontSize = 8;
        %hold on;
        % Extract unique InitCond values
        cond = subset.InitCond{1};

        bioreps = unique(subset.BioRep);
        shift = linspace(-0.05, 0.05, length(bioreps));
        for j=1:length(bioreps)
                biorep = bioreps(j);
                ssubset = subset(subset.BioRep==biorep,:);
                % Extract VariableValue, V0, and V0_err
                x = ssubset.VariableValue/1e3;
                %disp('X:');
                %disp(x);
                y = ssubset.V0;
                %disp('Y:');
                %disp(y);
                y_err = ssubset.V0_err;
                %disp('Y err:');
                %disp(y_err);
                errorbar(x+shift(j), y, y_err, aes{j}, 'DisplayName', sprintf('Init.cond., mkM: %s', cond));
                %errorbar(x, y, y_err, aes{j}, 'DisplayName', sprintf('Init.cond., mkM: %s', cond));
                hold on;
                %if (sign(y(end)<0))
                %        ylim([y(end)-0.5, 0]);
                %else
                %        ylim([0, y(end)+0.5]);
                %end
                %xlim([0 (max(x)+1)]);
                %disp(max(x)+1);
                xlabel([char(var), ', mM'], 'Interpreter', 'none', 'FontSize', fs);
                ylabel('V0, mkM/sec', 'Interpreter', 'none', 'FontSize', fs);
                legend('Location', 'northoutside', 'Interpreter', 'none');
        end
        hold off;
        if pageIs1st
                %disp(plot_path);
                exportgraphics(fig, plot_path, 'BackgroundColor', 'none', 'ContentType', 'vector');
                pageIs1st = false;
        else
                exportgraphics(fig, plot_path, 'Append', true, 'BackgroundColor', ...
                                'none', 'ContentType', 'vector');
        end 
        close(fig);
        %hold off;
    end
    %exportgraphics(fig, [plot_path filesep 'test_plot.pdf'], 'BackgroundColor', 'none', 'ContentType', 'vector');
    close all;
end

function tbl = getKineticsData(file_path)
        
        opts = detectImportOptions(file_path, ...
                                   'FileType', 'text', ...
                                   'TextType', 'string', ...
                                   'Delimiter', ',', ...
                                   'ReadVariableNames', true, 'ReadRowNames', false);
        opts.VariableTypes = [repmat({'string'},[1 3]), repmat({'double'},[1 4])];
        tbl = readtable(file_path, opts);
end