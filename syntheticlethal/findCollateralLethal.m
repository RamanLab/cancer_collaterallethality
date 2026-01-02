function [Collateral_Lethal,DLpairs_with_CollateralLethal]=findCollateralLethal(DLNormal_Fishers_FilePath,SLCancer_Fishers_FilePath,tissueName)
%DLNormal_Fishers_FilePath = '/home/maziya/CollateralLethality-2022new/BronchusLung_DLNormal_FishersResults.csv';
%SLCancer_Fishers_FilePath='/home/maziya/CollateralLethality-2022new/BronchusLung_SLCancer_FishersResults.csv';

DLNormal_Fishers = readcell(DLNormal_Fishers_FilePath);
SLCancer_Fishers = readcell(SLCancer_Fishers_FilePath);

DLNormal = split(DLNormal_Fishers(2:end,1),';');
SLCancer = SLCancer_Fishers(2:end,1);

Rxn1={};Rxn2={};
for index = 1:length(SLCancer)
    rowsMatch1 = strcmp(SLCancer{index, 1}, DLNormal(:, 1));
    Rxn1 = [Rxn1; DLNormal(rowsMatch1, :)];
    rowsMatch2 = strcmp(SLCancer{index, 1}, DLNormal(:, 2));
    Rxn2 = [Rxn2; DLNormal(rowsMatch2, :)];
end
DLpairs_with_CollateralLethal = cat(1, Rxn1, Rxn2);
[uniqueElements, ~, idxDuplicates] = unique(DLpairs_with_CollateralLethal, 'stable');
Collateral_Lethal = uniqueElements(histcounts(idxDuplicates) > 1);
CL = table('Size',[size(DLpairs_with_CollateralLethal,1),3],'VariableTypes',["string","string","string"],'VariableNames',["CollateralLethal","DoubleLethalR1","DoubleLethalR2"]);
Collateral_LethalPadded = [Collateral_Lethal; num2cell(NaN(height(CL) - length(Collateral_Lethal), 1))];
Collateral_LethalPadded(cellfun(@(x) isnumeric(x) && isnan(x), Collateral_LethalPadded)) = {'NaN'};
CL.CollateralLethal(:,1) = Collateral_LethalPadded;
CL(:,2:end) = DLpairs_with_CollateralLethal;
writetable(CL,[tissueName, '_CollateralLethal.csv']);
end

