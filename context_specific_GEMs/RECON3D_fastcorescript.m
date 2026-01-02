% clear;
%initCobraToolbox;
load data;
expData = log2(expMatrix+2);
load('Recon3d.mat');
model=Recon3D;
gene_expr = zeros(2302,length(expData(1,:)));
display(['Retrieving metabolic gene expression levels']);

% Get gene expression levels
for f=1:1:length(expData(1,:))
    genes_numeric = cellfun(@double, genes, 'UniformOutput', true);
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
            if max_gene_expr_row(i) - min_gene_expr_row(i) > 7
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

model = changeRxnBounds(model,{'EX_glc_D[e]'},-2,'l');
model = changeRxnBounds(model,{'EX_glc_D[e]'},-0.5,'u');

% Initialize tissueModel_rxns as a cell array before the loop
tissueModel_rxns = cell(length(sample_id(:,1)), 1);
rxn_mat = zeros(length(model.rxns),length(sample_id(:,1)));

% Initialize waitbar
f = waitbar(0, 'Model building in progress...');
pw = PoolWaitbar(length(sample_id(:,1)), 'TCGA model building in progress...');

for p = 1:length(sample_id(:,1))
    reaction_levels(:,p) = gene_to_reaction_levels_local(model, gene_id, row_gene_expr(:,p), @min, @max);
    
    HCR_t = prctile(nonzeros(reaction_levels(:,p)),70);
    ZCR_t = find(reaction_levels(:,p) == 0);
    
    if sum(strcmp(sample_type(p), 'Normal')) > 0
        HCR = [model.rxns(reaction_levels(:, p) > HCR_t); {'biomass_reaction'}; {'DM_atp_c_'}];
    else
        HCR = [model.rxns(reaction_levels(:, p) > HCR_t); {'biomass_reaction'}; {'DM_atp_c_'}; {'LDH_L'}; {'EX_lac_L[e]'}];
    end
    
    ZCR = model.rxns(ZCR_t);
    HMCR = model.rxns(reaction_levels(:, p) < HCR_t & reaction_levels(:, p) > prctile(nonzeros(reaction_levels(:, p)), 40));
    LMCR = model.rxns(reaction_levels(:, p) < prctile(nonzeros(reaction_levels(:, p)), 40) & reaction_levels(:, p) > 0);
    
    costbase = ones(length(model.rxns),1);
    costbase(findRxnIDs(model,HMCR)) = 0.1;
    costbase(findRxnIDs(model,LMCR)) = 100;
    costbase(findRxnIDs(model,ZCR)) = 1000000;
    costbase(findRxnIDs(model,HCR)) = 0;
    
% Get tissue-specific model
tissueModel = fastCoreWeighted(findRxnIDs(model,HCR), model, costbase, 1e-6);
% Store the reactions directly
tissueModel_rxns{p} = HCR;  % Store the core reactions
% Update reaction matrix using findRxnIDs
rxnIndices = findRxnIDs(model, HCR);
rxn_mat(rxnIndices, p) = 1;

% Create currentModel by removing unused reactions
currentModel = model;
currentModel.rxns = model.rxns(rxnIndices);
currentModel.S = model.S(:, rxnIndices);
currentModel.lb = model.lb(rxnIndices);
currentModel.ub = model.ub(rxnIndices);
currentModel.c = model.c(rxnIndices);
if isfield(model, 'subSystems')
    currentModel.subSystems = model.subSystems(rxnIndices);
end
if isfield(model, 'rxnNames')
    currentModel.rxnNames = model.rxnNames(rxnIndices);
end
    % Save the individual model
    % Get the sample ID and sanitize it for use as a filename
    current_sample_id = sample_id{p};
    % Remove any invalid characters from the filename
    current_sample_id = regexprep(current_sample_id, '[^\w\-]', '_');
    
    % Create the full path for saving
    save_path = fullfile('tissue_models', [current_sample_id '.mat']);
    
    % Save the model with its sample ID
    model_data = struct();
    model_data.model = currentModel;
    model_data.sample_id = sample_id{p};
    model_data.reaction_levels = reaction_levels(:,p);
    save(save_path, '-struct', 'model_data');
    
    waitbar(p / length(sample_id(:, 1)), f);
    increment(pw);
end

close(f);
delete(pw);

% Initialize rxn_mat
rxn_mat = zeros(length(model.rxns), length(expMatrix(1,:)));

% Populate rxn_mat
for p = 1:length(expMatrix(1,:))
    % Ensure we're working with cell arrays of character vectors
    if ~isempty(tissueModel_rxns{p})
        rxn_indices = ismember(model.rxns, tissueModel_rxns{p});
        rxn_mat(rxn_indices, p) = 1;
    end
end

% Save a summary of all models
save('tissue_models/model_summary.mat', 'rxn_mat', 'sample_id', 'tissueModel_rxns');

clear tissueModel_rxns;

Jac = zeros(length(expMatrix(1,:)),length(expMatrix(1,:)));
pw = PoolWaitbar(length(expMatrix(1,:)),'Jaccard distance being estimated...');
L = size(Jac,2);

for i = 1:length(expMatrix(1,:))
    expMatrix1 = expMatrix;
    rxn_mat1 = rxn_mat;
    temp_model1 = removeRxns(model,setdiff(model.rxns,model.rxns(find(rxn_mat1(:,i)))));
    for j = 1:L
        rxn_mat2 = rxn_mat1;
        temp_model2 = removeRxns(model,setdiff(model.rxns,model.rxns(find(rxn_mat2(:,j)))));
        Jac(i,j) = modelJaccardIndex(temp_model1, temp_model2);
    end
    increment(pw)
end
delete(pw);
function reaction_levels = gene_to_reaction_levels_local( model, genes, levels, f_and, f_or )
% Convert gene expression levels to reaction levels using GPR associations.
% Level is NaN if there is no GPR for the reaction or no measured genes.
%
% INPUTS
% model - cobra model
% genes - gene names
% levels - gene expression levels
% f_and - function to replace AND
% f_or - function to replace OR
%
% OUTPUTS
% reaction_levels - reaction expression levels
%
% Author: Daniel Machado, 2013
 reaction_levels = zeros(length(model.rxns), 1);
 for i = 1:length(model.rxns)
 level = eval_gpr(model.grRules{i}, genes, levels, f_and, f_or);
 reaction_levels(i) = level;
 end
end
function [result, status] = eval_gpr(rule, genes, levels, f_and, f_or)
% Evaluate the expression level for a single reaction using the GPRs.
% Note: Computes the expression level even if there are missing measured
% values for the given rule. This implementation is a modified version of
% an implementation provided in [Lee et al, BMC Sys Biol, 2012]
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
function [gene_id, gene_expr] = findUsedGenesLevels_local(model,num)
 % Find gene expression levels from the measured data (mRNA-Seq)
 % Input:
 % model - COBRA model struct (trimmed gene suffices, strings)
 % num - Nx2 matrix where the first column are the Entrez gene IDs
 % (integers) and the second column are the corresponding 
 % expression levels (from the mRNA-seq experiment).
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
 % In case multiple expression levels for genes are found in
 % the data, add values (happens only rarely)
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

function [tissue_model, rescue, HCtoMC, HCtoNC, MCtoNC] = CORDA_mod(model,metTests,...
 HCR,HMCR,LMCR,ZCR,PRtoNP,constraint,constrainby,om,ntimes,nl)
changeCobraSolver('gurobi7');
h = waitbar(0,'initializing waitbar');
OT = model.rxns(~ismember(model.rxns,cat(1,HCR,HMCR,LMCR,ZCR)));
if (nargin < 7) || isempty(PRtoNP)
 PRtoNP = 2;
end
if (nargin < 8) || isempty(constraint)
 constraint = 1;
end
if (nargin < 9) || isempty(constrainby)
 constrainby = 'val';
end
if (nargin < 10) || isempty(om)
 om = 1e+04;
end
if (nargin < 11) || isempty(ntimes)
 ntimes = 5;
end
if (nargin < 12) || isempty(nl)
 nl = 1e-02;
end
%% Step 1
costbase = zeros(length(model.rxns),1);
PR = [HMCR; LMCR];
NP = ZCR;
ES = HCR;
PRids = findRxnIDs(model,PR);
NPids = findRxnIDs(model,NP);
PRpres = false(length(PR),1);
NPpres = false(length(NP),1);
% costbase(PRids) = sqrt(om);
costbase(findRxnIDs(model,HMCR)) = nthroot(om,4)*10;
costbase(findRxnIDs(model,LMCR)) = nthroot(om,2)*10;
costbase(NPids) = om;
EStodelf = false(length(ES),1);
EStodelb = false(length(ES),1);
HCtoMC = zeros(length(ES),length(PR));
HCtoNC = zeros(length(ES),length(NP));
for j = 1:length(ES)
 waitbar(j/length(ES),h,'Step 1 - finding support reactions')
 
 %If reaction is sink then add it
 rem = false;
 if findRxnIDs(model,ES{j}) == 0
 id = find(strcmp(metTests(:,1),ES{j}));
 model = addReaction(model,metTests{id,1},metTests{id,2});
 costbase = [costbase; 0];
 rem = true;
 end
 
 %Test see which reactions are needed forward
 model = changeObjective(model,ES{j},1);
 flux = corsoFBA2(model,'max',constraint,constrainby,costbase+(nl*floor(10000*rand(length(costbase),1))/10000));
 if abs(flux.f) > 1e-6
 PRpres(abs(flux.x(PRids)) > 1e-6) = true;
 NPpres(abs(flux.x(NPids)) > 1e-6) = true;
 HCtoMC(j,abs(flux.x(PRids)) > 1e-6) = 1;
 HCtoNC(j,abs(flux.x(NPids)) > 1e-6) = 1;
 for k = 1:(ntimes-1)
 flux = corsoFBA2(model,'max',constraint,constrainby,costbase+(nl*floor(10000*rand(length(costbase),1))/10000));
 PRpres(abs(flux.x(PRids)) > 1e-6) = true;
 NPpres(abs(flux.x(NPids)) > 1e-6) = true;
 HCtoMC(j,abs(flux.x(PRids)) > 1e-6) = 1;
 HCtoNC(j,abs(flux.x(NPids)) > 1e-6) = 1;
 end
 else
 EStodelf(j) = true;
 end
 
 %Test see which reactions are needed backwards
 if model.lb(findRxnIDs(model,ES{j})) < 0
 flux = corsoFBA2(model,'min',-constraint,constrainby,costbase+(nl*floor(10000*rand(length(costbase),1))/10000));
 if abs(flux.f) > 1e-6
 PRpres(abs(flux.x(PRids)) > 1e-6) = true;
 NPpres(abs(flux.x(NPids)) > 1e-6) = true;
 HCtoMC(j,abs(flux.x(PRids)) > 1e-6) = 1;
 HCtoNC(j,abs(flux.x(NPids)) > 1e-6) = 1;
 for k = 1:(ntimes-1)
 flux = corsoFBA2(model,'min',-constraint,constrainby,costbase+(nl*floor(10000*rand(length(costbase),1))/10000));
 PRpres(abs(flux.x(PRids)) > 1e-6) = true;
 NPpres(abs(flux.x(NPids)) > 1e-6) = true;
 HCtoMC(j,abs(flux.x(PRids)) > 1e-6) = 1;
 HCtoNC(j,abs(flux.x(NPids)) > 1e-6) = 1;
 end
 else
 EStodelb(j) = true;
 end
 else
 EStodelb(j) = true;
 end
 
 if EStodelf(j) && EStodelb(j)
 fprintf(['Reaction ' ES{j} ' was deleted\n'])
 end
 
 %Remove sink if sink was added
 if rem
 model = removeRxns(model,ES{j});
 costbase(end) = [];
 end
end
%Tailor model
EStodel = EStodelf & EStodelb;
fprintf([num2str(length(find(EStodel))) ' blocked reactions removed from ES\n'])
ES(EStodel) = [];
HCtoMC(EStodel,:) = [];
HCtoNC(EStodel,:) = [];
HCtoMC = mat2dataset(HCtoMC,'varnames',PR,'obsnames',ES);
HCtoNC = mat2dataset(HCtoNC,'varnames',NP,'obsnames',ES);
fprintf([num2str(length(find(PRpres))) ' PR reactions added to ES\n'])
ES = cat(1,ES,PR(PRpres));
PR(PRpres) = [];
fprintf([num2str(length(find(NPpres))) ' NP reactions added to ES\n'])
ES = cat(1,ES,NP(NPpres));
NP(NPpres) = [];
%% Step 2
%% Step 2.1
%Assign costs
costbase = zeros(length(model.rxns),1);
costbase(findRxnIDs(model,NP)) = om;
%Initialize variables
NPid = findRxnIDs(model,NP);
PRxNP = zeros(length(PR),length(NP));
PRtodelf = false(length(PR),1);
PRtodelb = false(length(PR),1);
for j = 1:length(PR)
 waitbar(j/length(PR),h,'Step 2.1 - Checking PR and NP co-occurence')
 model = changeObjective(model,PR{j},1);
 
 %Check forward flux
 flux = corsoFBA2(model,'max',constraint,constrainby,costbase+(nl*floor(10000*rand(length(costbase),1))/10000));
 if isempty(flux.f)
 PRtodelf(j) = true;
 else
 for k = 1:(ntimes-1)
 PRxNP(j,abs(flux.x(NPid)) > 1e-6) = 1;
 flux = corsoFBA2(model,'max',constraint,constrainby,costbase+(nl*floor(10000*rand(length(costbase),1))/10000));
 end
 PRxNP(j,abs(flux.x(NPid)) > 1e-6) = 1;
 end
 
 %Check backwards flux
 if model.lb(findRxnIDs(model,PR{j})) < 0
 flux = corsoFBA2(model,'min',-constraint,constrainby,costbase+(nl*floor(10000*rand(length(costbase),1))/10000));
 if isempty(flux.f)
 PRtodelb(j) = true;
 else
 for k = 1:(ntimes-1)
 PRxNP(j,abs(flux.x(NPid)) > 1e-6) = 1;
 flux = corsoFBA2(model,'min',-constraint,constrainby,costbase+(nl*floor(10000*rand(length(costbase),1))/10000));
 end
 PRxNP(j,abs(flux.x(NPid)) > 1e-6) = 1;
 end
 else
 PRtodelb(j) = true;
 end
 if PRtodelf(j) && PRtodelb(j)
 fprintf(['Reaction ' PR{j} ' was deleted.\n'])
 end
end
PRxNP(PRtodelf & PRtodelb,:) = [];
PR(PRtodelf & PRtodelb) = [];
MCtoNC = mat2dataset(PRxNP,'varnames',NP,'obsnames',PR);
%% Step 2.2
%Add high occuring NPs to PR
t = sum(PRxNP);
ind = NP(t >= PRtoNP);
fprintf([num2str(length(ind)) ' reactions from NP are added to ES\n'])
%Fix occurence matrix
PR = cat(1,PR,ind);
PRxNP = [PRxNP;zeros(length(ind),length(NP))];
PRxNP(:,ismember(NP,ind)) = [];
NP(ismember(NP,ind)) = [];
%See which reactions from PR are no longer feasible
model = changeRxnBounds(model,NP,0,'b');
PRtodelf = false(length(PR),1);
PRtodelb = false(length(PR),1);
res1 = {};
res2 = {};
for j = 1:length(PR)
 waitbar(j/length(PR),h,'Step 2.2 - Checking PR feasibility')
 model = changeObjective(model,PR{j},1);
 
 %Check forward flux
 flux = optimizeCbModel(model,'max');
 if abs(flux.f) < 1e-6
 PRtodelf(j) = true;
 end
 
 %Check backwards flux
 if model.lb(findRxnIDs(model,PR{j})) < 0
 flux = optimizeCbModel(model,'min');
 if abs(flux.f) < 1e-6
 PRtodelb(j) = true;
 end
 else
 PRtodelb(j) = true;
 end
 
 if PRtodelf(j) && PRtodelb(j)
 fprintf([PR{j} ' was deleted. Dependent on: '])
 res1 = cat(1,res1,PR{j});
 tmp = find(PRxNP(j,:));
 if isempty(tmp)
 % NP reaction added can be dependent on other NP reactions that
 % appear less than PRtoNP times, and have thus been blocked.
 fprintf('Undefined')
 res2 = cat(1,res2,' ');
 else
 for k = 1:length(tmp)
 fprintf(NP{tmp(k)})
 if k ~= length(tmp)
 fprintf(', ')
 end
 if k == 1
 res2 = cat(1,res2,NP{tmp(k)});
 else
 res2{length(res2)} = strcat(res2{length(res2)},',',NP{tmp(k)});
 end
 end
 end
 fprintf('\n')
 end
end
fprintf([num2str(length(find(PRtodelf & PRtodelb))) ' reactions deleted from PR\n'])
PR(PRtodelf & PRtodelb) = [];
ES = cat(1,ES,PR);
rescue = cat(2,res1,res2);
%% Step 3
%Block reactions not in ES or OT
model = changeRxnBounds(model,...
 model.rxns(~ismember(model.rxns,cat(1,ES,OT))),0,'b');
%Define cost
costbase = zeros(length(model.rxns),1);
costbase(findRxnIDs(model,OT)) = om;
%Parse through ES
OTid = findRxnIDs(model,OT);
ESxOT = zeros(length(ES),length(OT));
for j = 1:length(ES)
 waitbar(j/length(ES),h,'Step 3 - Define remaining reactions')
 %Add sink if sink is needed
 rem = false;
 if findRxnIDs(model,ES{j}) == 0
 id = find(strcmp(metTests(:,1),ES{j}));
 model = addReaction(model,metTests{id,1},metTests{id,2});
 costbase = [costbase; 0];
 rem = true;
 end
 model = changeObjective(model,ES{j},1);
 
 %optimize ES forward
 for k = 1:ntimes
 flux = corsoFBA2(model,'max',constraint,constrainby,costbase+((nl*floor(10000*rand(length(costbase),1))/10000)));
 ESxOT(j,abs(flux.x(OTid)) > 1e-6) = 1;
 end
 
 %optimize ES backwards
 if model.lb(findRxnIDs(model,ES{j})) < 0
 for k = 1:ntimes
 flux = corsoFBA2(model,'min',-constraint,constrainby,costbase+((nl*floor(10000*rand(length(costbase),1))/10000)));
 ESxOT(j,abs(flux.x(OTid)) > 1e-6) = 1;
 end
 end
 
 %Remove sink if sink was added
 if rem
 model = removeRxns(model,ES{j});
 costbase(end) = [];
 end
end
fprintf([num2str(length(find(sum(ESxOT)))) ' reactions added to ES for final model\n'])
ES = cat(1,ES,OT(sum(ESxOT) ~= 0));
tissue_model = removeRxns(model,model.rxns(~ismember(model.rxns,ES)));
close(h)
end
function flux = corsoFBA2(model,onstr,constraint,constrainby,costas)
if strcmp(constrainby,'perc')
 constraint = abs(constraint);
end
%determine objective function
flux1 = optimizeCbModel(model,onstr);
if abs(flux1.f) < 1e-6
% warning('FBA problem infeasible')
 flux.f = [];
 flux.x = zeros(length(model.rxns),1);
 return
end
%relax results to avoid computational error
if strcmp(constrainby,'perc')
 flux1.f = flux1.x(model.c ~= 0)*(constraint/100); %Bound
elseif strcmp(constrainby,'val')
 if (flux1.f < constraint) && strcmp(onstr,'max')
 error('Objective Flux not attainable')
 elseif (flux1.f > constraint) && strcmp(onstr,'min')
 error('Objective Flux not attainable')
 else
 flux1.f = constraint; %Bound
 end
else
 error('Invalid Constraint option');
end
%save original model
model1 = model;
%See if cost is of right length
if length(costas) == 1
 costas = ones(length(model.rxns),1);
end
if ~iscolumn(costas)
 costas = costas';
end
if length(costas)==length(model.rxns)
 costas = [costas; costas];
elseif length(costas) ~= 2*length(model.rxns)
 fprintf('Invalid length of costs\n');
 flux = [];
 return
end
%find internal reactions that are actively reversible
orlen = length(model.rxns);
leng = find(model.lb<0 & model.ub>=0); 
%Tailor model
model.S = [model.S -model.S(:,leng); sparse(zeros(1,orlen+length(leng)))];
model.mets{length(model.mets)+1} = 'pseudomet';
model.b = zeros(length(model.mets),1);
model.c = zeros(orlen+length(leng),1);
model.ub = [model.ub; -model.lb(leng)];
model.lb = zeros(orlen+length(leng),1);
model.S(end,:) = [costas(1:orlen); costas(orlen+leng)];
model.rxns = cat(1,model.rxns,strcat(model.rxns(leng),'added'));
%add reaction for pseudomet consumption
model.rxns = cat(1,model.rxns,'EX_pseudomet');
model.ub = [model.ub; 1e20];
model.lb = [model.lb; 0];
temp = zeros(length(model.mets),1);
temp(end) = -1;
model.S = [model.S temp];
model.c = [model.c; 1];
%change bounds on original optimized reaction
t = find(model1.c);
for k = 1:length(t)
 model = changeRxnBounds(model,model1.rxns(t(k)),flux1.f(k),'b');
 if findRxnIDs(model,[model1.rxns{t(k)} 'added']) ~= 0
 model = changeRxnBounds(model,[model1.rxns{t(k)} 'added'],...
 0,'b');
 end
end
%perform FBA
flux2 = optimizeCbModel(model,'min');
flux.x = flux2.x(1:orlen);
flux.x(leng) = flux.x(leng) - flux2.x((orlen+1):(end-1));
flux.x(abs(flux.x) < 1e-8) = 0;
if isfield(flux1,'y')
 flux.y = flux1.y;
end
if isfield(flux1,'f')
 flux.f = flux1.f;
end
if isfield(flux2,'f')
 flux.fm = flux2.f;
end
end
function J = modelJaccardIndex(model1, model2)
% Compute intersection
JI = intersect(model1.rxns,model2.rxns);
% Compute union
JU = union(model1.rxns, model2.rxns);
% Calculate Jaccard indices
J = (length(JI)/ length(JU));
end