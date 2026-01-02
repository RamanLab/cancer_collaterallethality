
csvFilePath = '/Collateral/BronchusLung_CollateralLethal_new.csv';  
modelsFolderPath = '/Collateral/BLT';  
tissueName = 'Lung';
data = readtable(csvFilePath);


reactionsToCheck = data.Delete;  
collateralGenes = data.CollateralLethal_symbol;  
collateralReactions = data.CollateralLethal;  
deleteSymbols = data.Delete_symbol;  

%Filter rows: both 'CollateralLethal_symbol' and 'Delete_symbol' should be non-empty
validIdx = ~(cellfun(@isempty, collateralGenes) | cellfun(@isempty, deleteSymbols));
validReactions = reactionsToCheck(validIdx);
validCollateralGenes = collateralGenes(validIdx);
validCollateralReactions = collateralReactions(validIdx);
validDeleteSymbols = deleteSymbols(validIdx);


modelFiles = dir(fullfile(modelsFolderPath, '*.mat'));
numModels = length(modelFiles);

%Initialize a map to count successful checks
successCounts = containers.Map(validReactions, zeros(size(validReactions)));


for i = 1:numModels
    modelFile = modelFiles(i).name;
    modelPath = fullfile(modelsFolderPath, modelFile);
    loadedModel = load(modelPath);
    
    if isfield(loadedModel, 'model')
        model = loadedModel.model;
        modelReactions = model.rxns;      
        
        for j = 1:length(validReactions)
            deleteReaction = validReactions{j};          
            collateralReaction = validCollateralReactions{j};    
            % Check: deleteReaction absent AND collateralReaction present
            if ~ismember(deleteReaction, modelReactions) && ismember(collateralReaction, modelReactions)
                successCounts(deleteReaction) = successCounts(deleteReaction) + 1;
            end
        end
    end
end

%Find reactions satisfying the criteria in >=80% models
threshold = 0.8 * numModels;  
absentReactions = {};
presentCollateralGenes = {};
absentDeleteSymbols = {};

for i = 1:length(validReactions)
    reaction = validReactions{i};
    if successCounts(reaction) >= threshold
        absentReactions{end+1,1} = reaction; 
        presentCollateralGenes{end+1,1} = validCollateralGenes{i};  
        absentDeleteSymbols{end+1,1} = validDeleteSymbols{i};  
    end
end

%Create the results table
numAbsent = length(absentReactions);
tissueColumn = repmat({tissueName}, numAbsent, 1);

resultsTable = table(tissueColumn, presentCollateralGenes, absentDeleteSymbols, ...
    'VariableNames', {'Tissue', 'CollateralLethal_symbol', 'Delete_symbol'});
resultsTable = unique(resultsTable);

%Save
outputFileName = sprintf('AbsentReactionsSummary_%s.csv', tissueName);
writetable(resultsTable, outputFileName);
