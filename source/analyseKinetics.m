
function [kavg] = analysekinetics(nva)

        arguments
                nva.kinetics_data_path char
                nva.output_dir         char
        end

        kinetics = getKineticsData(nva);
        % average tech-reps
        groups = findgroups(kinetics.ID);
        
        kavg = table();
        kavg.ID = unique(kinetics.ID);
        kavg.InitCond = splitapply(@(x) x(1), kinetics.InitCond, groups);
        kavg.Variable = splitapply(@(x) x(1), kinetics.Variable, groups);
        kavg.VariableValue = splitapply(@(x) x(1), kinetics.VariableValue, groups);
        kavg.V0 = splitapply(@mean, kinetics.V0, groups);
        kavg.V0_err = splitapply(@(x) (sqrt(sum(x.^2)))/length(x), kinetics.V0_err, groups);
        kavg.BioRep = splitapply(@(x) x(1), kinetics.BioRep, groups);

        InitCond = [];
        for i=1:height(kavg)
                InitCond = [InitCond; ...
                            regexprep(regexprep(kavg.InitCond(i), join([":?" kavg.Variable(i) "=[0-9]+:?"],""),":"), "(^:|:&)", "")];
        end
        kavg.InitCond = InitCond;

        writetable(kavg, [nva.output_dir filesep 'kinetics_avg.csv'], ...
                        'FileType', 'text', ...
                        'Delimiter', ',');
        
        
end

function tbl = getKineticsData(nva)
        
        opts = detectImportOptions(nva.kinetics_data_path, ...
                                   'FileType', 'text', ...
                                   'TextType', 'string', ...
                                   'Delimiter', ',', ...
                                   'ReadVariableNames', true, 'ReadRowNames', false);
        opts.VariableTypes = [repmat({'string'},[1 3]), repmat({'double'},[1 5])];
        tbl = readtable(nva.kinetics_data_path, opts);
end