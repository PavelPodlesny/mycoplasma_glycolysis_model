function SBtabToSBMLmodel(varargin)
%% SBtabToSBMLmodel(directory, model_name, save_as_matfile)
%% Summary. This function reads an SBtab document with model's data, converts it to a structure, buildes a model, verifies the model and writes .txt file with model's equations.
%% Input arguments.
%        directory       -- char    -- an directory with model's data
%        model_name      -- char    -- model's name
%        save_as_matfile -- logical -- save the structure as an matfile or not


    dummy=0;
    arguments = {'test', 'test_model', true};
    if ~isempty(varargin)
        arguments(1:nargin) = varargin;
    end
    
    directory = num2str(arguments{1});
    model_name = arguments{2};
    save_as_matfile = arguments{3};

    data_path = ['/mnt/ecell_data/data2/input' filesep directory];
    
    model_data = SBtabToStruct(directory, save_as_matfile);
    [model, equations] = buildSBMLmodel(model_data, model_name);

    %% save model's equations
    fid = fopen([data_path filesep 'equations.txt'], 'w');
    fprintf(fid, '%s', equations);
    fclose(fid);
    %% save the model
    sbmlexport(model, [data_path filesep model_name '.xml']);
    
end