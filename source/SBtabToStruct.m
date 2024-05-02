function model_data = SBtabToStruct(varargin)
%% model_data = SBtabToStruct(directory, save_as_matfile)
%% Summary. This function reads an SBtab document with model's data and converts it to a structure.
%% Input arguments.
% directory       -- char    -- an directory with model's data 
% save_as_matfile -- logical -- save the structure as an matfile or not
%% Output arguments.
% model_data      -- struct  -- the structure with model's data

    dummy=0;
    arguments = {'test', false};
    
    if ~isempty(varargin)
        arguments(1:nargin) = varargin;
    end
    
    directory = num2str(arguments{1});
    save_as_matfile = arguments{2};
    
    data_path = ['/mnt/ecell_data/data2/input' filesep directory];
    
    %% read SBtab tables
    tables_list = dir([data_path filesep '*.tsv']);
    model_data  = struct;
    for i=1:length(tables_list)
        table_name = tables_list(i).name;
        field_name = lower(table_name(1:end-4));
        opts = detectImportOptions([data_path filesep table_name], 'FileType', 'text', ...
                                   'TextType', 'string', 'Delimiter', '\t', ...
                                   'VariableNamingRule', 'preserve', 'NumHeaderLines', 1, ...
                                   'ReadVariableNames', true, 'ReadRowNames', false);
        table = readtable([data_path filesep table_name], opts);
        if strcmp(field_name,'reaction')
            for i=1:height(table)
                IsReversible = table{i, '!IsReversible'};
                substitution = "->";
                if IsReversible
                    substitution = "<->";
                end
                table{i, '!ReactionFormula'} = replace(table{i, '!ReactionFormula'}, "<=>", substitution);
            end
        end
        model_data.(field_name) = table;
    end
    
    if save_as_matfile
        save([data_path filesep 'model_data.mat'], '-v7.3', '-struct','model_data');
    end
end