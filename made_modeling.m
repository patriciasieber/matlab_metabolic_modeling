%% run analysis with MADE in combination with iMAT

%% author: Patricia Sieber

%% import cobra
cd('/opt/gurobi752/linux64/matlab/');
gurobi_setup
cd('/home/patricia/bin/cobratoolbox/');
initCobraToolbox

out_dir = '/data/';
cd(out_dir);    

addpath(genpath('scripts/')); % genpath for recursive recognition in directory
addpath(genpath('/home/patricia/bin/made110209/')); % genpath for recursive recognition in directory

cmpi.init() %% initialize the MILP solver
%cmpi.set_solver('gurobi') % use the Gurobi solver

%recon3d_path = '/data/pidomics/data/annotations/Recon3D/Recon3D_301.mat';
%load(recon3d_path);


atlas_path = '/data/annotations/atlas/humanGEM.mat';
atlas_out_dir = '/data/annotations/atlas/';

load(atlas_path);


%% load atlas annotation and load features into files

% Write gene list to file
fileID = fopen(strcat(atlas_out_dir, 'atlas', '_gene_ids.txt'), 'w');
for n = 1:length(ihuman.genes)
   fprintf(fileID, '%s\n', ihuman.genes{n});
end
fclose(fileID);


% Write metabolite names to file
fileID = fopen(strcat(atlas_out_dir, 'atlas', '_met_names.txt'), 'w');
for n = 1:length(ihuman.metNames)
    fprintf(fileID, '%s\n', ihuman.metNames{n});
end
fclose(fileID);


% Write stoichiometric matrix to file
S = full(ihuman.S);
fileID = fopen(strcat(atlas_out_dir, 'atlas', '_stoich_mat.txt'), 'w');
formatspec = [repmat('%.0f\t', 1, size(S, 2) - 1), '%.0f\n'];
for n = 1:size(S,1)
    fprintf(fileID, formatspec, S(n ,:));
end
fclose(fileID);


% Write reaction rules and properties to file
fileID = fopen(strcat(atlas_out_dir, 'atlas', '_gr_rules.txt'), 'w');
fprintf(fileID, '%s\n', 'Gene reaction rule');
for n = 1:length(ihuman.grRules)
    fprintf(fileID, '%s\n', ihuman.grRules{n});
end
fclose(fileID);

% tested also for Recon3D slim, and Recon2.2 but does not work (due to entrez ids)


%% run MADE
m = make_elf_model(ihuman);

% compare IA vs control
% use genes from DEG analysis
% FC of MRN values, p-value (FDR) from DESeq2

patients = {'P01','P02','P03','P04','P05','P06'};
gene_names = importdata('/data/made/genenames_onlyatlas.txt');

timepoints={'early','middle','late'};

for i= 1:length(patients)
       
    fold_change = table2array(readtable(strcat('/data/made/patient_',patients{i},'_tps_foldchanges_onlyatlas.csv')));
    p = table2array(readtable(strcat('/data/made/patient_',patients{i},'_tps_pvalues_onlyatlas.csv')));
    
    % run MADE
    [gene_states,genes,sol,models] = made(m,fold_change,p,1/3, ...
                                        'gene_names',gene_names)

     save(strcat('/data/made/',patients{i},'results_atlas_made.mat'),'gene_states','genes','sol','models');
     
     % write into separate files as input for imat
     mrn = table2array(readtable(strcat('/data/made/patient_',patients{i},'_tps_mrn_onlyatlas.csv')));
     
     
     %for j= 1:size(mrn,2)
     for j=1:length(timepoints)
        % gene id, expression value, minus, plus
        
        minus = zeros(length(gene_names),1);   
       
        imat_input = table(gene_names, mrn(:,j), minus, gene_states(:,j));
        imat_input.Properties.VariableNames = {'Gene_id';'Values';'minus';'plus'}; 
     
        %writetable(imat_input, strcat('/data/made/preprocessing/patient_', patients{i}, '_tp', num2str(j), '_preprocessed_made.txt'), 'Delimiter','\t');
        writetable(imat_input, strcat('/data/made/preprocessing/patient_', patients{i}, '_' , timepoints{j}, '_preprocessed_made.txt'), 'Delimiter','\t');

        
        imat_input = table(gene_names, gene_states(:,j));
        imat_input.Properties.VariableNames = {'Gene_id';'Values'};
        %writetable(imat_input, strcat('/data/made/preprocessing/patient_', patients{i}, '_tp', num2str(j), '_preprocessed_made_onlyValues.txt'), 'Delimiter','\t');
        writetable(imat_input, strcat('/data/made/preprocessing/patient_', patients{i}, '_' , timepoints{j}, '_preprocessed_made_onlyValues.txt'), 'Delimiter','\t');

     end
     
                              
end

% does not run easily for Recon annotations (due to entrez ids) 
% continue with human metabolic atlas id

%% run imat
changeCobraSolver('gurobi','ALL');

% Load base human metabolism model Recon3D
model = ihuman;  % does not work with original ihuman annotation but with its elf-model m
recon_genes = model.genes;

data_dir = '/data/made/preprocessing/';
out_dir = 'results/models_iMAT_made/';
% Make sure the output directory exists
if (exist(out_dir, 'dir') == 0)
    system(char(strcat('mkdir -p',{' '},out_dir)));  %% {' '} for empty space
end

% Get data files
file_list = strsplit(ls(data_dir));
fl_take = find(contains(file_list, 'made'));
file_list = file_list(fl_take);
file_out_list = strrep(file_list, 'txt', 'mat');

% Filter input files processed in a prior run
file_junk = cellfun(@(x) exist(x, 'file') == 2, strcat(out_dir, file_out_list));
file_out_list(file_junk) = [];
file_list(file_junk) = [];

% Parse GPRs - use later for mapping of gene expression to reaction expression
%model = initFBAFields(model);  % rules were missing in ihuman model

model = convertOldStyleModel(model);  % to insert rules, needed for next steps
gprs = GPRparser(model); 


parfor (n = 1:length(file_list), 2)  %% 2 for number of parallel processes
    expr_data = struct(); % need to initialize structs for parallel optimization
    options = struct();
    gdat = struct();
    %[expr_data.Locus, expr_data.Data, expr_data.minus, expr_data.plus] = import_expression_file(strcat(data_dir,file_list{n}));
    [expr_data.Locus, expr_data.Data] = import_expression_file(strcat(data_dir,file_list{n}));
    
    % Convert the locus field to cell if necessary
    if (~iscell(expr_data.Locus))
        expr_data.Locus = cellstr(expr_data.Locus);
    end
    changeCobraSolver('gurobi','ALL'); % necessary for parallel optimization
    [gpres, gloc] = ismember(recon_genes, expr_data.Locus);
    strcat(num2str(sum(gpres)), ' model genes matched to expression gene set')
    expr_data.Data = expr_data.Data + abs(min(expr_data.Data)); % avoids propagation of negative expression values from Microarray data to reaction expression, where COBRA expects values == -1 or >= 0
    gdat.gene = model.genes(gpres);
    gdat.value = expr_data.Data(gloc(gpres));
    options.solver = 'iMAT';
    options.expressionRxns = my_mapExpressionToReaction(model, gprs, gdat); 
    if (length(model.rxns) ~= length(options.expressionRxns))
        error('Malformed reaction expression vector; quitting');
    end
    if (sum(options.expressionRxns > -1) == 0)
        error('Undefined expression in every reaction');
    end
    %options.threshold_lb = max(expr_data.Data(expr_data.minus == 1));
        
    %options.threshold_lb = -1;   % in 1sttrial with 0, 2nd trial with -1
    %options.threshold_ub = min(expr_data.Data(expr_data.plus == 1));
    
    % Try with o and 1 as values
    options.threshold_lb = 0;
    options.threshold_ub = 1;
    
    % Create the costum model
    % excluded modelCom from iMAT.m file in /home/patricia/bin/cobratoolbox/src/dataIntegration/transcriptomics/iMAT/iMAT.m
    tissueModel = createTissueSpecificModel(model, options);
    % Save costum model
    mf = matfile(strcat(out_dir, file_out_list{n}), 'writable', true); %for usage in parfor
    mf.Rxns = tissueModel.rxns;
    %Rxns = tissueModel.rxns;
    %save(strcat(out_dir, file_out_list{n}), 'Rxns');
end

%% iMAT reaction extraction
% Get unique list of reaction ids
rxn_ref_list = model.rxns;

% Get list of subsystem classification
subsys_class = model.subSystems;
if (iscell(subsys_class{1}))
    subsys_class = cellfun(@(x) x{1}, subsys_class, 'UniformOutput', false);
end

% Select directories for algorithm
data_imat_dir = 'results/models_iMAT_made/';
out_dir_imat = 'results/reactionData_iMAT_made/';

% Make sure the output directory exists
if (exist(out_dir_imat, 'dir') == 0)
    mkdir(out_dir_imat);
end

% Get data files
%file_list = strsplit(ls(data_imat_dir), {'\n', '\t'});
file_list = strsplit(ls(data_imat_dir), {'\n', ' '});
fl_take = contains(file_list, '.mat');
file_list = file_list(fl_take);
file_out_list = strrep(file_list, 'mat', 'csv');

% Filter input files processed in a prior run
file_junk = ismember(file_out_list, strsplit(ls(out_dir_imat), {'\n', '\t'}));
file_out_list(file_junk) = [];
file_list(file_junk) = [];

% Write reaction list to disk along with subsystem mebmership
fileID = fopen(strcat(out_dir_imat, 'rxn_list.csv'),'wt');
formatSpec = '%s\t%s\n';
fprintf(fileID,formatSpec,'Rxn_id,', 'Subsystem');
[nrows] = length(rxn_ref_list);
for row = 1:nrows
    fprintf(fileID,formatSpec, rxn_ref_list{row}, subsys_class{row});
end
fclose(fileID);
clear fileID formatSpec nrows row

% Extract reaction information
% Rxns.* corresponds to tissueModel.rxns
parfor (n = 1:length(file_list),3)
    mf = matfile(strcat(data_imat_dir, file_list{n}));
        if (isstruct(mf.Rxns))
            Rxns = mf.Rxns;
        else
            Rxns = struct();
            Rxns.Expressed = mf.Rxns;
            Rxns.DownRegulated = {};
            Rxns.UpRegulated = {};
            Rxns.UnExpressed = {};
            Rxns.UnknownIncluded = {};
        end
    % Set operations on reference reacton list
    rxn_matrix = zeros(length(rxn_ref_list), 5);
    rxn_matrix(:,1) = ismember(rxn_ref_list, Rxns.DownRegulated);
    rxn_matrix(:,2) = ismember(rxn_ref_list, Rxns.UpRegulated);
    rxn_matrix(:,3) = ismember(rxn_ref_list, Rxns.UnknownIncluded);
    rxn_matrix(:,4) = ismember(rxn_ref_list, Rxns.Expressed);
    rxn_matrix(:,5) = ismember(rxn_ref_list, Rxns.UnExpressed);
    % Save as csv
    fileID = fopen(strcat(out_dir_imat, file_out_list{n}),'wt');
    formatSpec = '%s\t%s\t%s\t%s\t%s\n';
    fprintf(fileID,formatSpec, 'DownRegulated', 'UpRegulated', 'UnknownIncluded', 'Expressed', 'UnExpressed');
    formatSpec = '%u\t%u\t%u\t%u\t%u\n';
    [nrows,~] = size(rxn_matrix);
    fprintf(fileID,formatSpec,rxn_matrix');
    fclose(fileID);
end

%% iMAT reaction analysis
% with R
%system('R CMD BATCH scripts/analyseRxns_edit.R iMAT'); %% TODO: change name of input file in R script to start it from here!

