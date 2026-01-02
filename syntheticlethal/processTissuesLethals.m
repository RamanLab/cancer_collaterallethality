function [LethalCounts] = processTissuesLethals(tissue_identifier_normal, tissue_identifier_cancer, lethalrxnlist, lethaltype)
% tissue_identifier_normal = 'BLN[1-9]';
% tissue_identifier_cancer = 'BLT[1-9]';
%lethalrxnlist is a cell array where each column contains the single lethal reactions from one
%model. This output is obtained from fastSL.
%lethaltype denotes 'single' or 'double' lethals being given in the
%lethalrxnlist
% Process normal tissues
[lethals_normal,number_normal] = processTissue(tissue_identifier_normal, lethalrxnlist,lethaltype);
counts_normal = countAndCreateArray(lethals_normal);

% Process cancer tissues
[lethals_cancer,number_cancer] = processTissue(tissue_identifier_cancer, lethalrxnlist,lethaltype);
counts_cancer = countAndCreateArray(lethals_cancer);

% Combine unique lethal values from normal and cancer tissues
LethalCounts = unique([counts_normal(:, 1); counts_cancer(:, 1)]);

% Update LethalCounts array for normal and cancer tissues
LethalCounts = updateLethalCounts(LethalCounts, counts_normal, lethals_normal, 2);
LethalCounts = updateLethalCounts(LethalCounts, counts_cancer, lethals_cancer, 3);

tf = cellfun('isempty', LethalCounts);
LethalCounts(tf) = {0};
switch lethaltype
    case 'single'
        LethalCounts = [LethalCounts, num2cell([number_normal, number_cancer] - (cell2mat(LethalCounts(:,2:3))))];
        LethalCounts = cell2table(LethalCounts);
        LethalCounts.Properties.VariableNames = {'SingleLethal','Counts in Normal','Counts in Cancer','Not SL in Normal','Not SL in Cancer'};
        writetable(LethalCounts, ['SL_LethalCounts_', regexp(tissue_identifier_normal, '.*(?=N)', 'match', 'once'), '.csv']);
    case 'double'
        rowsToDelete = cellfun(@(x) strcmp(x, ';'), LethalCounts(:, 1));
        LethalCounts(rowsToDelete, :) = [];
        LethalCounts = [LethalCounts, num2cell([number_normal, number_cancer] - (cell2mat(LethalCounts(:,2:3))))];
        LethalCounts = cell2table(LethalCounts);
        LethalCounts.Properties.VariableNames = {'DoubleLethal','Counts in Normal','Counts in Cancer','Not SL in Normal','Not SL in Cancer'};
        writetable(LethalCounts, ['DL_LethalCounts_', regexp(tissue_identifier_normal, '.*(?=N)', 'match', 'once'), '.csv']);
end
end

function [lethals, numberofmodels] = processTissue(tissue_identifier, lethalrxnlist,lethaltype)
matches = cellfun(@(x) ~isempty(regexp(x, tissue_identifier, 'all')), lethalrxnlist(1, :))';
indices = find(matches);
switch lethaltype
        case 'single'
            lethals ={};
            for str = 1:length(indices)
                lethals = [lethals, lethalrxnlist(2:end, indices(str))];       
            end
            numberofmodels = numel(indices);
        case 'double'
            indices = indices(1:2:end-1,:);
            lethals = cell(length(lethalrxnlist)-1,length(indices));
            index =1;  
            for str = 1:length(indices)
                lethals(:,index) = lethalrxnlist(2:end, indices(str));
                lethals(:,index+1) =lethalrxnlist(2:end, indices(str)+1);
                lethals(:,index)  = strcat(lethals(:,index),';',lethals(:,index+1));
                lethals(:,index+1) =[];
                index = index+1;
            end
             numberofmodels = numel(indices);
end
end

function counts = countAndCreateArray(lethals)
lethals = lethals(~cellfun('isempty', lethals));
UniqLethalRxn = unique(lethals);
lethal_counts = cellfun(@(x) sum(ismember(lethals, x)), UniqLethalRxn, 'UniformOutput', false);
counts = [UniqLethalRxn, lethal_counts];
end

function LethalCounts = updateLethalCounts(LethalCounts, counts, lethals, columnIndex)
[common_elements, idx] = ismember(LethalCounts(:, 1), counts(:, 1));
valid_indices = common_elements & idx ~= 0;
LethalCounts(valid_indices, columnIndex) = counts(idx(valid_indices), 2);
end