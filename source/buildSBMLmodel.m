function [model, equations] = buildSBMLmodel(model_data, varargin)
%% sbml_model = buildSBMLmodel(model_data, model_name)
%% Summary. This function builds SBML model from 'model_data' structure and verifies it.
%% Input arguments.
% model_data -- struct  -- an structure containing model's data
% model_name -- char    -- model's name
%% Output arguments.
% model     -- SimBiology model object -- a model
% equations -- char -- model's equations 
    dummy=0;
    arguments = {'mycoplasma_glycolysis'};
    
    if ~isempty(varargin)
        arguments(1:nargin-1) = varargin;
    end

    model_name = arguments{1};
    
    %%create model object 
    model = sbiomodel(model_name);

    %% add compartments
    for i=1:height(model_data.compartment)
        comp = addcompartment(model, model_data.compartment{i, '!Name'});
        comp.Capacity = model_data.compartment{i, '!Size'};
        comp.CapacityUnits = model_data.compartment{i, '!Unit'};
    end
    %% add parameters to the model (parameters are shared across all model's entities)
    for i=1:height(model_data.parameter)
        param = addparameter(model, model_data.parameter{i, '!ID'});
        param.Value = model_data.parameter{i, '!DefaultValue'};
        param.Units = model_data.parameter{i, '!Unit'};
        param.Notes = model_data.parameter{i, '!Name'};
        param.Tag   = model_data.parameter{i, '!QuantityType'};
    end
    %% add species
    for i=1:height(model_data.compound)
        comp = sbioselect(model, 'Type', 'compartment', 'Name', model_data.compound{i, '!Location'});
        species = addspecies(comp, model_data.compound{i, '!ID'});
        species.Value = model_data.compound{i, '!InitialValue'};
        species.Units = model_data.compound{i, '!Unit'};
        species.Notes = model_data.compound{i, '!Name'};
    end
    %% add reaction
    for i=1:height(model_data.reaction)
        rct =  addreaction(model, model_data.reaction{i, '!ReactionFormula'});
        rct.ReactionRate = model_data.reaction{i, '!KineticLaw'};
        rct.Reversible = model_data.reaction{i, '!IsReversible'};
        rct.Notes = model_data.reaction{i, '!Name'};
        rct.Tag   = model_data.reaction{i, '!KineticLaw:Name'};
    end
    %% verify model
    verify(model);
    %% get model's equations
    equations = getequations(model);
end