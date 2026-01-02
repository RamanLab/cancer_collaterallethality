%Identifying active and inactive reactions after removing the reaction from the normal model in 
%order to examine whether the resulting flux distribution resembles that of the cancer model
inputFolderPath = '/Collateral/FluxDiffResults_LIN';  
outputFolderPath = '/Collateral/FluxDiffResults_LIN/onlyactiveinactive';  
if ~exist(outputFolderPath, 'dir')
    mkdir(outputFolderPath);
end
csvFiles = dir(fullfile(inputFolderPath, '*.csv'));
for i = 1:length(csvFiles)
    filePath = fullfile(inputFolderPath, csvFiles(i).name);
    data = readtable(filePath);
    data = rmmissing(data);
    filteredData = data((data.FluxCancer == 0 & data.FluxNormal ~= 0) | (data.FluxNormal == 0 & data.FluxCancer ~=0), :);
    [~, fileName, ext] = fileparts(csvFiles(i).name);
    outputFileName = fullfile(outputFolderPath, [fileName, '_filtered', ext]);
    writetable(filteredData, outputFileName);
    fprintf('Processed and saved: %s\n', outputFileName);
end
