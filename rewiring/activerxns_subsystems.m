
%listing out all the rxns that are on or off in cancer and normal-del along
%with subsystems
model = readCbModel('Recon.mat');  
inputFolder = '/Collateral/FluxDiffResults_PN/onlyactiveinactive';  
csvFiles = dir(fullfile(inputFolder, '*.csv'));
allData = table();

for i = 1:length(csvFiles)
    filePath = fullfile(inputFolder, csvFiles(i).name);
    data = readtable(filePath);
    allData = [allData; data];
end

uniqueReactions = unique(allData.Reaction);

onInCancer = {};
offInCancer = {};
onInNormal = {};
offInNormal = {};
subsystemsCancerOn = {};
subsystemsCancerOff = {};
subsystemsNormalOn = {};
subsystemsNormalOff = {};


for i = 1:length(uniqueReactions)
    rxn = uniqueReactions{i}; 
    idx = strcmp(allData.Reaction, rxn);    
    fluxCancer = allData.FluxCancer(idx);
    fluxNormal = allData.FluxNormal(idx);
    
    % Determine ON/OFF status in Cancer
    if any(fluxCancer ~= 0)
        onInCancer{end+1,1} = rxn;
        rxnIdx = find(strcmp(model.rxns, rxn), 1);
        if ~isempty(rxnIdx)
            subsystemsCancerOn{end+1,1} = model.subSystems{rxnIdx};
        else
            subsystemsCancerOn{end+1,1} = 'Unknown';
        end
    else
        offInCancer{end+1,1} = rxn;
        rxnIdx = find(strcmp(model.rxns, rxn), 1);
        if ~isempty(rxnIdx)
            subsystemsCancerOff{end+1,1} = model.subSystems{rxnIdx};
        else
            subsystemsCancerOff{end+1,1} = 'Unknown';
        end
    end
    
    % Determine ON/OFF status in Normal
    if any(fluxNormal ~= 0)
        onInNormal{end+1,1} = rxn;
        rxnIdx = find(strcmp(model.rxns, rxn), 1);
        if ~isempty(rxnIdx)
            subsystemsNormalOn{end+1,1} = model.subSystems{rxnIdx};
        else
            subsystemsNormalOn{end+1,1} = 'Unknown';
        end
    else
        offInNormal{end+1,1} = rxn;
        rxnIdx = find(strcmp(model.rxns, rxn), 1);
        if ~isempty(rxnIdx)
            subsystemsNormalOff{end+1,1} = model.subSystems{rxnIdx};
        else
            subsystemsNormalOff{end+1,1} = 'Unknown';
        end
    end
end


onCancerTable = table(onInCancer, subsystemsCancerOn, 'VariableNames', {'Reaction', 'Subsystem'});
offCancerTable = table(offInCancer, subsystemsCancerOff, 'VariableNames', {'Reaction', 'Subsystem'});
onNormalTable = table(onInNormal, subsystemsNormalOn, 'VariableNames', {'Reaction', 'Subsystem'});
offNormalTable = table(offInNormal, subsystemsNormalOff, 'VariableNames', {'Reaction', 'Subsystem'});

outputFolder = fullfile(inputFolder, 'fluxdiff_rxns_subs');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%Write to CSV files
writetable(onCancerTable, fullfile(outputFolder, 'onInCancer.csv'));
writetable(offCancerTable, fullfile(outputFolder, 'offInCancer.csv'));
writetable(onNormalTable, fullfile(outputFolder, 'onInNormal.csv'));
writetable(offNormalTable, fullfile(outputFolder, 'offInNormal.csv'));

