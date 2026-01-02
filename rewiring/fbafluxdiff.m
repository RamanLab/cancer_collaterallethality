filePath = '/CL_Validationfba.csv'; 
data = readtable(filePath, 'ReadVariableNames', true);
reactionsToDelete = data.('RxnToDelete');
cancerModels = data.('cancerModel');
normalModels = data.('normalModel');
cancerTypes = data.('Cancer');  

basePath = '/Collateral';


uniqueCancers = unique(cancerTypes);

% Create an output folder for each unique Cancer
outputFolders = containers.Map;
for i = 1:numel(uniqueCancers)
    cancerType = uniqueCancers{i};
    outputFolder = fullfile(basePath, sprintf('FluxDiffResults_%s', cancerType));
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    outputFolders(cancerType) = outputFolder;
end


for i = 1:height(data)
    cancerName = cancerModels{i};
    normalName = normalModels{i};
    reactionToDelete = reactionsToDelete{i};
    cancerType = cancerTypes{i}; 

    outputFolder = outputFolders(cancerType);

    cancerFolder = findExactModelFolder(basePath, cancerName);
    if isempty(cancerFolder)
        warning('Cancer model folder not found for %s', cancerName);
        continue;
    end
    cancerModelPath = fullfile(basePath, cancerFolder, [cancerName, '.mat']);
    if exist(cancerModelPath, 'file')
        cancerModel = load(cancerModelPath);
        cancerModel = cancerModel.model;    
    else
        warning('Cancer model file not found: %s', cancerModelPath);
        continue;
    end


    normalFolder = findExactModelFolder(basePath, normalName);
    if isempty(normalFolder)
        warning('Normal model folder not found for %s', normalName);
        continue;
    end
    normalModelPath = fullfile(basePath, normalFolder, [normalName, '.mat']);
    if exist(normalModelPath, 'file')
        normalModel = load(normalModelPath);
        normalModel = normalModel.model;
    else
        warning('Normal model file not found: %s', normalModelPath);
        continue;
    end

    % Delete reaction in normal model
    normalModel = changeRxnBounds(normalModel, reactionToDelete, 0, 'b');
    solCancer = optimizeCbModel(cancerModel);
    solNormal = optimizeCbModel(normalModel);
    
    solutionCancer = solCancer.x;
    solutionNormal = solNormal.x;
    
    allReactions = union(cancerModel.rxns, normalModel.rxns, 'stable');
    fluxCancer = nan(length(allReactions), 1);
    fluxNormal = nan(length(allReactions), 1);

    [~, idxCancer] = ismember(cancerModel.rxns, allReactions);
    fluxCancer(idxCancer) = solutionCancer;
    [~, idxNormal] = ismember(normalModel.rxns, allReactions);
    fluxNormal(idxNormal) = solutionNormal;
    fluxDifference = fluxCancer - fluxNormal;
    fluxDifference(isnan(fluxCancer) | isnan(fluxNormal)) = NaN;


    outputTable = table(repmat({cancerName}, length(allReactions), 1), ...
                        repmat({normalName}, length(allReactions), 1), ...
                        repmat({reactionToDelete}, length(allReactions), 1), ...
                        allReactions, fluxCancer, fluxNormal, fluxDifference, ...
                        'VariableNames', {'CancerModel', 'NormalModel', 'DeletedReaction', 'Reaction', 'FluxCancer', 'FluxNormal', 'FluxDifference'});


    outputFileName = sprintf('FluxDiff_results_%s_vs_%s.csv', cancerName, normalName);
    outputFilePath = fullfile(outputFolder, outputFileName);
    writetable(outputTable, outputFilePath);

    fprintf('Processed %s (%s) vs %s | Deleted: %s | Saved: %s\n', cancerName, cancerType, normalName, reactionToDelete, outputFilePath);
end
