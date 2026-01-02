% clear;
%initCobraToolbox;
load data;
expData = log2(expMatrix+2);
load('recon3Econsistent','model');
gene_expr = zeros(2697,length(expData(1,:)));
display(['Retrieving metabolic gene expression levels']);
for f=1:1:length(expData(1,:))
    [gene_id, gene_expr(:,f)] = findUsedGenesLevels_local(model, [genes_numeric, expData(:,f)]);

end
clear f;

display(['Preprocessing the expression data to determine expressed metabolic genes']);
thr = zeros(length(gene_id),1);
mean_gene_expr_row = zeros(length(gene_id),1);
pr_25_gene_expr_row = zeros(length(gene_id),1);
pr_90_gene_expr_row = zeros(length(gene_id),1);
min_gene_expr_row = zeros(length(gene_id),1);
max_gene_expr_row = zeros(length(gene_id),1);

for i=1:1:length(gene_id)
    mean_gene_expr_row(i) = mean(nonzeros(gene_expr(i,:)));
    pr_25_gene_expr_row(i) = prctile(nonzeros(gene_expr(i,:)),25);
    pr_90_gene_expr_row(i) = prctile(nonzeros(gene_expr(i,:)),90);
    if ~(isempty(min(nonzeros(gene_expr(i,:)))))
        min_gene_expr_row(i) = min(nonzeros(gene_expr(i,:)));
    end
    if ~(isempty(max(nonzeros(gene_expr(i,:)))))
        max_gene_expr_row(i) = max(nonzeros(gene_expr(i,:)));
    end
    mean_gene_expr_row(isnan(mean_gene_expr_row))=0;
    pr_25_gene_expr_row(isnan(pr_25_gene_expr_row))=0;
    pr_90_gene_expr_row(isnan(pr_90_gene_expr_row))=0;
    min_gene_expr_row(isnan(min_gene_expr_row))=0;
    max_gene_expr_row(isnan(max_gene_expr_row))=0;
    if pr_90_gene_expr_row(i) < 1.586
        thr(i) = 1.586;
    else
        if min_gene_expr_row(i) > 1.586
            if  max_gene_expr_row(i) - min_gene_expr_row(i) > 7
                thr(i) = mean_gene_expr_row(i);
            else
                thr(i) = min_gene_expr_row(i);
            end
        else
                thr(i) = mean_gene_expr_row(i);    
        end
    end
end
clear i max_gene_expr_row min_gene_expr_row mean_gene_expr_row pr_90_gene_expr_row pr_25_gene_expr_row;

row_gene_expr = zeros(length(gene_id),length(expData(1,:)));
for a=1:1:length(gene_id)
       for j=1:1:length(expData(1,:))
            if gene_expr(a,j) >= thr(a)
                row_gene_expr(a,j) = 5*log(1+(gene_expr(a,j)/thr(a)));
            else
                row_gene_expr(a,j) = 0;
            end
       end       
end
clear a j;

model = changeRxnBounds(model,{'EX_glc(e)'},-2,'l');
model = changeRxnBounds(model,{'EX_glc(e)'},-0.5,'u');

tissueModel_rxns = cell(length(sample_id(:,1)), 1);
rxn_mat = zeros(length(model.rxns),length(sample_id(:,1)));
tissue_models = cell(length(sample_id(:,1)), 1); % Store all tissue models

f = waitbar(0, 'Model building in progress...');
pw = PoolWaitbar(length(sample_id(:,1)), 'TCGA model building in progress...');

for p = 1:length(sample_id(:,1))
    reaction_levels(:,p) = gene_to_reaction_levels_local(model, gene_id, row_gene_expr(:,p), @min, @max);
    
    HCR_t = prctile(nonzeros(reaction_levels(:,p)),70);
    ZCR_t = find(reaction_levels(:,p) == 0);
    
    if sum(strcmp(sample_type(p), 'Normal')) > 0
        HCR = [model.rxns(reaction_levels(:, p) > HCR_t); {'biomass_reaction'}; {'DM_atp[c]'}];
    else
        HCR = [model.rxns(reaction_levels(:, p) > HCR_t); {'biomass_reaction'}; {'DM_atp[c]'}; {'LDH_L'}; {'EX_lac_L(e)'}];
    end
    
    ZCR = model.rxns(ZCR_t);
    HMCR = model.rxns(reaction_levels(:, p) < HCR_t & reaction_levels(:, p) > prctile(nonzeros(reaction_levels(:,p)), 40));
    LMCR = model.rxns(reaction_levels(:, p) < prctile(nonzeros(reaction_levels(:,p)), 40) & reaction_levels(:, p) > 0);
    
    weights = ones(length(model.rxns), 1);
    weights(findRxnIDs(model,HMCR)) = 0.1;
    weights(findRxnIDs(model,LMCR)) = 100;
    weights(findRxnIDs(model,ZCR)) = 1000000;
    weights(findRxnIDs(model,HCR)) = 0;
    
    coreInd = findRxnIDs(model, HCR);
    
    %[reconstruction, reconInd, LP] = swiftcore(model, coreInd, weights, 1e-5, true, 'gurobi');
    try
    [reconstruction, reconInd, LP] = swiftcore(model, coreInd, weights, 1e-5, true, 'gurobi');
    
    % Save results for debugging
    save(sprintf('debug_model_%d.mat', p), 'model', 'coreInd', 'weights', 'reconstruction', 'reconInd', 'LP');
    
    % Check if flux is empty or not assigned
    if isempty(LP)
        disp(['Flux is empty for model: ', num2str(p)]);
        
        % Save info in the workspace for later analysis
        assignin('base', sprintf('failed_model_%d', p), model);
        assignin('base', sprintf('failed_LP_%d', p), LP);
        assignin('base', sprintf('failed_coreInd_%d', p), coreInd);
        assignin('base', sprintf('failed_weights_%d', p), weights);
        
        continue;  % Skip this iteration if flux is empty
    end
    
catch ME
    warning('Error in model %d: %s', p, ME.message);
    
    % Save error info in the workspace for debugging
    assignin('base', sprintf('error_model_%d', p), model);
    assignin('base', sprintf('error_message_%d', p), ME.message);
    
    continue;  % Skip this iteration if there was an error in swiftcore
end

    % Inside the loop where you're processing the models
%try
    %[reconstruction, reconInd, LP] = swiftcore(model, coreInd, weights, 1e-6, true, 'gurobi');
    
    % Check if flux is empty or not assigned
    %if isempty(LP)  % You may need to check LP or flux depending on what core returns
        %disp('Flux is empty for model: ');
        %disp(model);  % Print model information if flux is empty
        %continue;  % Skip this iteration if flux is empty
    %end
    
%catch ME
    %warning('Error in model %d: %s', p, ME.message);
    %continue;  % Skip this iteration if there was an error in swiftcore
%end

    
    tissueModel_rxns{p} = reconstruction.rxns;
    rxn_mat(reconInd, p) = 1;
    tissue_models{p} = reconstruction; % Store the complete model
    
    waitbar(p / length(sample_id(:,1)), f);
    increment(pw);
end
close(f);
delete(pw);

% Calculate Jaccard indices
Jac = zeros(length(expMatrix(1,:)), length(expMatrix(1,:)));
pw = PoolWaitbar(length(expMatrix(1,:)), 'Jaccard distance being estimated...');

for i = 1:length(expMatrix(1,:))
    expMatrix1 = expMatrix;
    rxn_mat1 = rxn_mat;
    model1 = tissue_models{i};
    
    for j = 1:length(expMatrix(1,:))
        model2 = tissue_models{j};
        Jac(i,j) = modelJaccardIndex(model1, model2);
    end
    increment(pw);
end
delete(pw);

% Save all results
save('tissue_specific_models.mat', 'tissue_models', 'rxn_mat', 'Jac', 'tissueModel_rxns', 'reaction_levels');

% Helper functions remain the same
function [gene_id, gene_expr] = findUsedGenesLevels_local(model,num)
    genes_ID = zeros(length(model.genes),1);
    for i = 1:length(model.genes)
        genes_ID(i) = str2num(model.genes{i});
    end
    cnts = -1*ones(length(genes_ID),2);
    cnts(:,1) = genes_ID;
    for i = 1:length(genes_ID)
        cur_ID = genes_ID(i);
        flag = 0;
        for ii = 1:length(num)
            if num(ii,1) == cur_ID
                if flag == 1
                    cnts(i,2) = cnts(i,2) + num(ii,2);
                end
                if flag == 0
                    flag = 1;
                    cnts(i,2) = num(ii,2);
                end
            end
        end
    end
    data_inds = find(cnts(:,2)~= -1);
    gene_expr = cnts(data_inds,2);
    gene_id = model.genes(data_inds);
end

function reaction_levels = gene_to_reaction_levels_local(model, genes, levels, f_and, f_or)
    reaction_levels = zeros(length(model.rxns), 1);
    for i = 1:length(model.rxns)
        level = eval_gpr(model.grRules{i}, genes, levels, f_and, f_or);
        reaction_levels(i) = level;
    end
end

function [result, status] = eval_gpr(rule, genes, levels, f_and, f_or)
    EVAL_OK = 1;
    PARTIAL_MEASUREMENTS = 0;
    NO_GPR_ERROR = -1;
    NO_MEASUREMENTS = -2;
    MAX_EVALS_EXCEEDED = -3;
    MAX_EVALS = 1000;
    NONETYPE = 'NaN';
    NUMBER = '[0-9\.\-e]+';
    MAYBE_NUMBER = [NUMBER '|' NONETYPE];
    
    expression = rule;
    result = NaN;
    status = EVAL_OK;
    
    if isempty(expression)
        status = NO_GPR_ERROR;
    else
        rule_genes = setdiff(regexp(expression,'\<(\w|\-)+\>','match'), {'and', 'or'});
        total_measured = 0;
        
        for i = 1:length(rule_genes)
            j = find(strcmp(rule_genes{i}, genes));
            if isempty(j)
                level = NONETYPE;
            else
                level = num2str(levels(j));
                total_measured = total_measured + 1;
            end
            expression = regexprep(expression, ['\<', rule_genes{i}, '\>'], level );
        end
        
        if total_measured == 0
            status = NO_MEASUREMENTS;
        else
            if total_measured < length(rule_genes)
                status = PARTIAL_MEASUREMENTS;
            end
            
            maybe_and = @(a,b)maybe_functor(f_and, a, b);
            maybe_or = @(a,b)maybe_functor(f_or, a, b); 
            str_wrapper = @(f, a, b)num2str(f(str2double(a), str2double(b)));
            counter = 0;
            
            while isnan(result)
                counter = counter + 1;
                if counter > MAX_EVALS
                    status = MAX_EVALS_EXCEEDED;
                    break
                end
                try 
                    result = eval(expression);            
                catch e   
                    paren_expr = ['\(\s*(', MAYBE_NUMBER,')\s*\)'];
                    and_expr = ['(',MAYBE_NUMBER,')\s+and\s+(',MAYBE_NUMBER,')'];
                    or_expr = ['(',MAYBE_NUMBER,')\s+or\s+(',MAYBE_NUMBER,')'];
                    expression = regexprep(expression, paren_expr, '$1');
                    expression = regexprep(expression, and_expr, '${str_wrapper(maybe_and, $1, $2)}');
                    expression = regexprep(expression, or_expr, '${str_wrapper(maybe_or, $1, $2)}');
                end
            end
        end
    end
end

function c = maybe_functor(f, a, b)
    if isnan(a) && isnan(b)
        c = nan;
    elseif ~isnan(a) && isnan(b)
        c = a;
    elseif isnan(a) && ~isnan(b)
        c = b;
    else 
        c = f(a,b);
    end
end

function J = modelJaccardIndex(model1, model2)
    JI = intersect(model1.rxns, model2.rxns);
    JU = union(model1.rxns, model2.rxns);
    J = (length(JI) / length(JU));
end