function [Fisher_results] = doFishersLethal(LethalCountsFilePath,tissueName,lethaltype)
%Input contigency table of normal and cancer single and double lethals
%use Fishers test to identify SL significant in cancer and DL significant in normal models
%tissueName = 'Breast';
% LethalCountsFilePath = '/home/maziya/CollateralLethality-2022new/SL_LethalCounts_BL.csv';
LethalCountsTable = readtable(LethalCountsFilePath);
pvalue_threshold = 0.05;
switch lethaltype
    case 'single'
        for i = 1:size(LethalCountsTable,1)
            x= table([LethalCountsTable{i,2};LethalCountsTable{i,4}],[LethalCountsTable{i,3};LethalCountsTable{i,5}],'VariableNames',{'Normal','Cancer'},'RowNames',{'Lethal','NotLethal'});
            [h(i),p(i,1),stats(i)] = fishertest(x,'Tail','left','Alpha',0.01);
        end
        OddsRatio = [stats.OddsRatio].';
        ConfidenceInterval = {stats.ConfidenceInterval}.';
        ConfidenceInterval = cell2mat(ConfidenceInterval);
        Fisher_results = horzcat(LethalCountsTable,(array2table(horzcat(p,OddsRatio, ConfidenceInterval))));
        Fisher_results.Properties.VariableNames(6:9) ={'p-value','OddsRatio','ConfidenceInterval-I','ConfidenceInterval-II'};
        Fisher_results  = Fisher_results(Fisher_results.("p-value") <= pvalue_threshold, :);
        Fisher_results = sortrows(Fisher_results,6,'ascend');
        writetable(Fisher_results,[tissueName,'_SLCancer_FishersResults.csv']);
    case 'double'
        for i = 1:size(LethalCountsTable,1)
            x= table([LethalCountsTable{i,2};LethalCountsTable{i,4}],[LethalCountsTable{i,3};LethalCountsTable{i,5}],'VariableNames',{'Normal','Cancer'},'RowNames',{'Lethal','NotLethal'});
            [h(i),p(i,1),stats(i)] = fishertest(x,'Tail','right','Alpha',0.01);
        end
        OddsRatio = [stats.OddsRatio].';
        ConfidenceInterval = {stats.ConfidenceInterval}.';
        ConfidenceInterval = cell2mat(ConfidenceInterval);
        Fisher_results = horzcat(LethalCountsTable,(array2table(horzcat(p,OddsRatio, ConfidenceInterval))));
        Fisher_results.Properties.VariableNames(6:9) ={'p-value','OddsRatio','ConfidenceInterval-I','ConfidenceInterval-II'};        
        Fisher_results  = Fisher_results(Fisher_results.("p-value") <= pvalue_threshold, :);
        Fisher_results = sortrows(Fisher_results,6,'ascend');
        writetable(Fisher_results,[tissueName,'_DLNormal_FishersResults.csv']);
end
end