function SBtabToSBMLmodel(nva)
%% SBtabToSBMLmodel(directory, model_name, save_as_matfile)
%% Summary. This function reads an SBtab document with model's data, converts it to a structure, buildes a model, verifies the model and writes .txt file with model's equations.
%% Input arguments.
%        directory       -- char    -- a directory with model's data
%        model_name      -- char    -- model's name
%        save_as_matfile -- logical -- save the structure as an matfile or not

    arguments
            nva.directory                char = 'test'
            nva.model_name               char = 'test_model'
            nva.save_as_matfile       logical = true
            nva.unit_convers          logical = false
    end
    
    data_path = ['/mnt/ecell_data/data2/input' filesep nva.directory];
    
    model_data = SBtabToStruct(nva.directory, nva.save_as_matfile);
    [model, equations] = buildSBMLmodel(model_data, nva.model_name, nva.unit_convers);

    %% save model's equations
    fid = fopen([data_path filesep 'equations.txt'], 'w');
    fprintf(fid, '%s', equations);
    fclose(fid);
    %% save the model
    sbmlexport(model, [data_path filesep nva.model_name '.xml']);
    
end