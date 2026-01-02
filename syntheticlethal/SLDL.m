folder = 'models_wo_outliers';
model_list = get_model_names('models_wo_outliers');
model_list = model_list';

tic;
parfor n = 1:numel(model_list)
    %     changeCobraSolver('ibm_cplex','all');
    m = model_list{n};
    modelOrig = load(m);
    model = modelOrig.model;
    atpm = 'DM_atp[c]'; cutoff = 0.01;
    model.lb(find(model.c)) = 0.001;
    model = changeObjective(model,atpm);
    %     lac_id = findRxnIDs(model,'EX_lac_L(e)');
    %     if lac_id ~= 0
    %         model = changeObjective(model,'EX_lac_L(e)');
    %         soln = optimizeCbModel(model);
    %         if soln.stat ==1
    [selExc, ~] = findExcRxns(model);
    eliList= model.rxns(selExc==1);
    disp(n)
    try
        [Jsl{n,1},Jdl{n,1}]=parallelSL(model,cutoff,eliList);
    catch
        Jsl{n,1} ={};
        Jdl{n,1} ={};
    end
    %         else
    %             Jsl{n,1} ={};
    %             Jdl{n,1} ={};
    %         end
    %     end
    %   [sgd{n,1},dgd{n,1}]=parallelSL_Gene(model,cutoff,eliList);
end
toc;

% check row sizes
sizes_SL = cellfun('length', Jsl);
maxlen_SL = max(sizes_SL);
sizes_DL = cellfun('length', Jdl);
maxlen_DL = max(sizes_DL);

SL = cell2mat(cellfun(@(x) [x(:, 1); zeros(maxlen_SL - size(x, 1), 1)], Jsl, 'UniformOutput', false));
SL = reshape(SL(1:maxlen_SL * length(Jsl)), maxlen_SL, length(Jsl));

DL = zeros(maxlen_DL,length(Jdl));
column_no = 1;
for q = 1:length(Jdl)
    for r = 1:length(Jdl{q,1})
        DL(r,column_no:column_no+1) = Jdl{q,1}(r,:);
    end
    column_no = column_no+2;
end

SL_model = cell(size(SL));
for n = 1:numel(model_list)
    m = model_list{n};
    modelOrig = load(m);
    model = modelOrig.model;
    nonZeroIndices = SL(:, n) ~= 0;
    SL_model(nonZeroIndices, n) = model.rxns(SL(nonZeroIndices, n));
end

DL_model = cell(size(DL));
column_no = 1;
for n = 1:numel(model_list)
    m = model_list{n};
    modelOrig = load(m);
    model = modelOrig.model;
    nonZeroIndices = DL(:, column_no) ~= 0;
    DL_model(nonZeroIndices, column_no:column_no+1) = cellfun(@(x) model.rxns(x), num2cell(DL(nonZeroIndices, column_no:column_no+1)), 'UniformOutput', false);
    column_no = column_no + 2;
end

SL_model_new = cell(size(SL_model,1)+1,size(SL_model,2));
SL_model_new(2:end,1:end) = SL_model;
modelnames = model_list';
SL_model_new(1,1:end) = modelnames;
writecell(SL_model_new,'EXlac_L_SLmodels.csv')

DL_model_new = cell(size(DL_model,1)+1,size(DL_model,2));
DL_model_new(2:end,1:end) = DL_model;
modelname = repelem(model_list,2)';
DL_model_new(1,1:end) = modelname;
writecell(DL_model_new,'EXlac_L_DLmodels.csv')


%findinSLDLcountsg cumulative occurences of rxn in either normal or cancer models
% A is rxn and occurences in all normal or cancer models
% B is unique list of rxns
% rxnCount=[];
% for i = 1:length(B)
%     row = strcmp(B{i,1},A(:,1));
%     if ~isempty(row)
%         id = find(row==1);
%         for j = 1:length(id)
%             rxnCount(j,1) = (A{id(j),2});
%         end
%         rxnsum{i,1} = sum(rxnCount);
%     end
% end