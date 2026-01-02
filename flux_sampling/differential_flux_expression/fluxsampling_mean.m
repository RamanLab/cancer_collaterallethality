% Set the base parameters
folder = 'D:\Rutuja\FLUXSAMPLING_inputs';  % Base folder path
normal_identifier = 'N';  % Identifier for normal tissue models
cancer_identifier = 'T';  % Identifier for cancer tissue models

% Get list of tissue folders
tissue_folders = dir(folder);
tissue_folders = tissue_folders([tissue_folders.isdir]); % Get only directories
tissue_folders = tissue_folders(~ismember({tissue_folders.name}, {'.', '..'})); % Remove . and ..

% Process each tissue folder
for tissue_idx = 1:length(tissue_folders)
    tissue_name = tissue_folders(tissue_idx).name;
    tissue_path = fullfile(folder, tissue_name);
    
    disp(['Processing tissue: ' tissue_name]);
    
    try
        % Call the main function for each tissue
        [recon3e_data_rand_N, recon3e_data_rand_C, stats, pVals, rxn_meanflux_cancer, ...
            rxn_meanflux_normal, rxn_meanflux_cancer_samp, rxn_meanflux_normal_samp] = ...
            fluxsampling_auto(tissue_path, normal_identifier, cancer_identifier);
    catch ME
        % Print error information but continue processing
        disp(['Warning: Error in processing ' tissue_name]);
        disp(['Error message: ' ME.message]);
        disp('Continuing with next tissue...');
    end
end

function [recon3e_data_rand_N, recon3e_data_rand_C, stats, pVals, rxn_meanflux_cancer, ...
    rxn_meanflux_normal, rxn_meanflux_cancer_samp, rxn_meanflux_normal_samp] = ...
    fluxsampling_auto(folder, model_identifier_normal, model_identifier_cancer)
    
    % Load recon3E model
    disp('Loading recon3E model...');
    recon3e = load('recon3Econsistent_exch.mat');
    recon3e_rxns = recon3e.model.rxns;
    
    % Debug information
    disp(['Type of recon3e_rxns: ' class(recon3e_rxns)]);
    if iscell(recon3e_rxns)
        disp(['Size of recon3e_rxns: ' num2str(size(recon3e_rxns))]);
        disp(['Type of first element: ' class(recon3e_rxns{1})]);
    end

    % Get model names
    model_list = get_model_names(folder)';
    disp(['Found ' num2str(length(model_list)) ' model files']);

    % Count models for each condition
    count_models_N = sum(contains({model_list{:}}, model_identifier_normal));
    count_models_C = sum(contains({model_list{:}}, model_identifier_cancer));
    
    disp(['Normal models: ' num2str(count_models_N)]);
    disp(['Cancer models: ' num2str(count_models_C)]);

    % Initialize arrays
    recon3e_data_C = zeros(size(recon3e_rxns, 1), 2000 * min(count_models_C,count_models_N));
    recon3e_data_N = zeros(size(recon3e_rxns, 1), 2000 * min(count_models_C,count_models_N));

    % Initialize arrays for storing individual sample data
    normal_files = {};
    cancer_files = {};
    normal_data_per_sample = cell(1, count_models_N);
    cancer_data_per_sample = cell(1, count_models_C);
    normal_idx = 1;
    cancer_idx = 1;

    % Extract and store flux data
    count_N = 1;
    count_C = 1;
    for i = 1:length(model_list)
        disp(['Processing file ' num2str(i) ' of ' num2str(length(model_list))]);
        [~, filename, ~] = fileparts(model_list{i});
        
        opts = detectImportOptions(model_list{i});
        model_sample_org = readtable(model_list{i}, opts);
        model_sample_org(1,:) = [];
        model_sample_rxnlist = table2cell(model_sample_org(:, 1));
        
        % Debug information for model_sample_rxnlist
        disp(['Type of model_sample_rxnlist: ' class(model_sample_rxnlist)]);
        if ~isempty(model_sample_rxnlist)
            disp(['Type of first element in model_sample_rxnlist: ' class(model_sample_rxnlist{1})]);
        end
        
        is_normal_model = contains(model_list{i}, model_identifier_normal);
        
        try
            [~, idx_recon3e] = ismember(model_sample_rxnlist, recon3e_rxns);
        catch ME
            disp('Error in ismember operation:');
            disp(['Error message: ' ME.message]);
            disp('Attempting to continue...');
            continue;
        end
        
        if any(idx_recon3e)
            flux_data = zeros(length(recon3e_rxns), 2000);
            flux_data(idx_recon3e, :) = cell2mat(table2cell(model_sample_org(:, 2:end)));
            
            if is_normal_model
                recon3e_data_N(idx_recon3e, count_N:count_N+1999) = cell2mat(table2cell(model_sample_org(:, 2:end)));
                normal_files{normal_idx} = filename;
                normal_data_per_sample{normal_idx} = flux_data;
                normal_idx = normal_idx + 1;
                count_N = count_N + 2000;
            else
                recon3e_data_C(idx_recon3e, count_C:count_C+1999) = cell2mat(table2cell(model_sample_org(:, 2:end)));
                cancer_files{cancer_idx} = filename;
                cancer_data_per_sample{cancer_idx} = flux_data;
                cancer_idx = cancer_idx + 1;
                count_C = count_C + 2000;
            end
        end
    end

    % Random sampling (Original logic)
    disp('Performing random sampling...');
    samp = 2000;
    init = 2;
    for i = 1:min(count_models_C,count_models_N)
        rall{i, 1} = randi([init, samp], 100, 1);
        init = samp + 1;
        samp = samp + 1999;
    end
    rtotal = vertcat(rall{:});
    
    % Get random samples
    recon3e_data_rand_N = [];
    recon3e_data_rand_C = [];
    for index = 1:length(rtotal)
        recon3e_data_rand_N = [recon3e_data_rand_N recon3e_data_N(:, rtotal(index))];
        recon3e_data_rand_C = [recon3e_data_rand_C recon3e_data_C(:, rtotal(index))];
    end

    % Calculate means (Original logic)
    rxn_meanflux_cancer = mean(recon3e_data_rand_C, 2);
    rxn_meanflux_normal = mean(recon3e_data_rand_N, 2);
    rxn_meanflux_cancer_samp = mean(recon3e_data_C, 2);
    rxn_meanflux_normal_samp = mean(recon3e_data_N, 2);

    % Statistical analysis
    disp('Performing statistical analysis...');
    [stats, pVals] = compareTwoSamplesStat(recon3e_data_rand_C, recon3e_data_rand_N, ...
        {'ks','chiSquare','rankSum','tTest'});

    % Create results directory
    [~, tissue_name] = fileparts(folder);
    results_dir = fullfile(folder, 'results');
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end

    % Calculate individual sample means
    normal_sample_means = zeros(length(recon3e_rxns), length(normal_files));
    cancer_sample_means = zeros(length(recon3e_rxns), length(cancer_files));
    
    for i = 1:length(normal_files)
        normal_sample_means(:,i) = mean(normal_data_per_sample{i}, 2);
    end
    
    for i = 1:length(cancer_files)
        cancer_sample_means(:,i) = mean(cancer_data_per_sample{i}, 2);
    end

    % Save results
    disp('Saving results...');
    save_detailed_results(results_dir, recon3e_rxns, normal_sample_means, cancer_sample_means, ...
        rxn_meanflux_normal, rxn_meanflux_cancer, rxn_meanflux_normal_samp, rxn_meanflux_cancer_samp, ...
        normal_files, cancer_files, tissue_name);
    
    disp(['Completed processing for tissue: ' tissue_name]);
end

function save_detailed_results(results_dir, rxns, normal_sample_means, cancer_sample_means, ...
    normal_rand_means, cancer_rand_means, normal_all_means, cancer_all_means, ...
    normal_files, cancer_files, tissue_name)
    
    % Start with reactions column
    T = table(rxns, 'VariableNames', {'Reaction'});
    
    % Add columns for each normal sample
    for i = 1:length(normal_files)
        col_name = normal_files{i};
        T.(col_name) = normal_sample_means(:, i);
    end
    
    % Add columns for each cancer sample
    for i = 1:length(cancer_files)
        col_name = cancer_files{i};
        T.(col_name) = cancer_sample_means(:, i);
    end
    
    % Add overall means
    T.Normal_Random_Mean = normal_rand_means;  % From random sampling
    T.Cancer_Random_Mean = cancer_rand_means;  % From random sampling
    T.Normal_All_Mean = normal_all_means;      % From all samples
    T.Cancer_All_Mean = cancer_all_means;      % From all samples
    
    % Save to CSV
    filename = fullfile(results_dir, [tissue_name '_detailed_means.csv']);
    writetable(T, filename);
    disp(['Saved detailed means to: ' filename]);
end

function model_list = get_model_names(folder)
    % Get all CSV files in the folder (non-recursive)
    files = dir(fullfile(folder, '*.csv'));
    model_list = cellstr(fullfile({files.folder}, {files.name}));
end

function [stats, pVals] = compareTwoSamplesStat(data1, data2, testTypes)
    % Initialize output structures
    stats = struct();
    pVals = struct();
    
    % Ensure data is in correct format (columns are samples)
    if size(data1, 2) == 1
        data1 = data1';
    end
    if size(data2, 2) == 1
        data2 = data2';
    end
    
    % Perform each requested statistical test
    for i = 1:length(testTypes)
        switch lower(testTypes{i})
            case 'ks'
                % Kolmogorov-Smirnov test
                [h, p, ks_stat] = kstest2(data1', data2');
                stats.ks = ks_stat;
                pVals.ks = p;
                
            case 'chisquare'
                % Chi-square test
                [h, p, chi_stat] = chi2test(data1, data2);
                stats.chiSquare = chi_stat;
                pVals.chiSquare = p;
                
            case 'ranksum'
                % Wilcoxon rank sum test
                [p, h, stats_temp] = ranksum(data1', data2');
                stats.rankSum = stats_temp.ranksum;
                pVals.rankSum = p;
                
            case 'ttest'
                % Two-sample t-test
                [h, p, ~, stats_temp] = ttest2(data1', data2');
                stats.tTest = stats_temp.tstat;
                pVals.tTest = p;
        end
    end
end

function [h, p, chi_stat] = chi2test(data1, data2)
    % Simple implementation of chi-square test
    % Combine data and create frequency distributions
    binEdges = linspace(min([data1(:); data2(:)]), max([data1(:); data2(:)]), 50);
    
    counts1 = histcounts(data1(:), binEdges);
    counts2 = histcounts(data2(:), binEdges);
    
    % Calculate chi-square statistic
    expected = (counts1 + counts2) / 2;  % Expected counts under null hypothesis
    chi_stat = sum((counts1 - expected).^2 ./ expected) + ...
               sum((counts2 - expected).^2 ./ expected);
    
    % Calculate p-value (degrees of freedom = number of bins - 1)
    df = length(binEdges) - 2;
    p = 1 - chi2cdf(chi_stat, df);
    
    % Set significance level at 0.05
    h = (p < 0.05);
end