clear; clc;
% This analysis script can be used to run a parellelised version of 
% connectome-based predictive modelling with LOOCV or repeated k-fold CV.
% Authors: Rory Boyle & Yihe Weng
% Emails: rorytboyle@gmail.com & wengy@tcd.ie
% Date: 07/10/2021
% Date edited: 20/05/2022
%% 1) Load data and prepare variables
% Specify file paths
target_path = 'W:\CPM_FU3\ICV\doc\behav_data.csv';
covars_path = 'W:\CPM_FU3\ICV\doc\site_covariates.csv';
conn_mx_dir = 'W:\CPM_FU3\ICV\connectivity_matrices';

thresholds = [0.1,0.05,0.01,0.005,0.001,0.0005,0.0001];
All_pos = zeros(length(thresholds),1);
All_neg = zeros(length(thresholds),1);
All_combined = zeros(length(thresholds),1);

for t = 1: length(thresholds)
 
    output_path = ['W:\Yihe\CPM_FU3\ICV\output\',['Diffthrshold_both',num2str(thresholds(t))]];
    if ~exist(output_path,'dir')
        mkdir(output_path)
    end
    
    % Specify target variable (ensure string is same as listed in .csv file)
    target_var_name = 'ICV';
    
    % Specify name of model for saving output
    model_name = 'ICV_CPM';
    
    % Prepare input data
    nsubs = 1120;  % list number of participants
    [cpm_predictors, file_order] = prep_predictors_CPM(conn_mx_dir, nsubs,...
        output_path);
    
    % Load predictor variables (i.e. connectivity matrices)
    all_mats = cpm_predictors;
    
    % Prepare target variable 
    [cpm_target, final_ppts, ix_ppts_to_keep] = prep_target_CPM(target_path,...
        file_order, output_path);
    
    % Load target variable and covariates to be included
    data = readtable(covars_path);
    
    % Get ppts with connectivity matrices, covariates, and target variable 
    subids=table2array(data(:,1));
    log_ix = ismember(subids, final_ppts);
    data = data(log_ix,:);  % drop ppts w/o conn mx/covariates/target variable
    
    % Specify target variable
    all_behav = cpm_target;
    
    % Prep covariates - age, sex, fwd: covar_names = {'age', 'sex', 'mean_FWD'};
    % if no covars to be included, use: covar_names = {}
    % covariates must be named here as listed in .csv file
    covar_names = {'sex','age','meanFD','site1','site2','site3','site4','site5','site6','site7'};
    all_covars = table2array(data(:, covar_names));
    
    % convert covar_names to string for saving
    covar_str = "";
    for i = 1:length(covar_names)
        if i == 1
            covar_str = append(covar_str, covar_names(i));
        else
            covar_str = append(covar_str, ' ', covar_names(i));
        end
    end
    
    %% 2) Specify model inputs
    % Specify method for dealing with confounds. adjust_stage = 'relate' 
    % adjusts for confounds via partial correlation during feature selection.
    % adjust_stage = 'fit' adjusts for confounds via multiple regression during
    % model fitting. adjust_stage = 'both' adjusts for confounds at both above
    % steps. adjust_stage = '' does not adjust for confounds.
    adjust_stage = 'both';
    
    % Specify cross-validation scheme - LOOCV: k=n, 10-fold CV: k=10, etc.
    % for LOOCV: k = length(all_behav);
    k = 10;
    
    % Specify iterations (may want to run multiple iterations of k-fold CV)
    % for LOOCV: iterations = 1
    iterations=50;  
    
    % Specify feature selection type. 'p-value' will threshold edges based on
    % p-value of correlation between edges and target variable. 'sparsity' will
    % threshold edges by selecting the X % most strongly correlated edges to 
    % target variable.
    thresh_type = 'p-value';
    
    % Specify feature selection threshold. If 'p-value', enter p-value for
    % correlations between edges and target variable. If 'sparsity', enter % of
    % most highly correlated edges to be retained. 
    % Note: % should be in decimal (i.e. between 0 and 1, 5% = 0.05)
    thresh = thresholds(t);
    
    % Specify edge frequency threshold (i.e. how many folds edge must be
    % significantly correlated (i.e. below thresh) with target variable in
    % order to be selected for application to the test set)
    % Note: % should be in decimal (i.e. between 0 and 1, 100% = 1)
    freq_thresh = 1;
    
    % Specify if permutation test to be conducted (will greatly increase
    % runtime of code); run_permutation = 'yes' or 'no'
    run_permutation = 'no';
    
    % if yes, specify number of iterations of random permutation 
    % (e.g. perm_iterations = 1000;). if no, save as []
    perm_iterations = 1000;
    
    % Save info describing model inputs 
    model_info = struct('target_path',target_path,'covars_path',covars_path, ...
        'target_variable', target_var_name, 'N', length(all_behav),'covars',...
        covar_str, 'adjust_stage', {adjust_stage}', 'k', k, 'iterations',...
        iterations, 'fs_thresh_type',thresh_type,'fs_thresh',thresh,'freq_thresh',...
        freq_thresh,'permutation_test', run_permutation, 'perm_iterations',...
         perm_iterations, 'output_path', output_path, 'model_name', model_name);
    
    clearvars -except all_behav all_mats all_covars k thresh thresh_type freq_thresh...
        adjust_stage iterations model_info perm_iterations run_permutation thresholds ...
    target_path covars_path conn_mx_dir t All_pos All_neg All_combined
    
    %% 3) Preallocate arrays for storing results and parameters 
    % Preallocate arrays for storing CPM predicted values and results
    [behav_pred_pos_all, behav_pred_neg_all, behav_pred_combined_all, ...
        R_pos, P_pos, R_neg, P_neg, R_combined, P_combined,...
        rsq_pos, rsq_neg, rsq_combined, mae_pos, mae_neg, mae_combined] = ...
        prep_results_arrays_CPM(all_mats, iterations);
    
    % Preallocate arrays for storing CPM parameters
    if strcmp(adjust_stage, 'relate')
        [~, ~, int_pos_ntwrk, int_neg_ntwrk,...
        int_combined_ntwrk, slope_pos_ntwrk, slope_neg_ntwrk, ...
        slope_combined_ntwrk, ~, ~, ~] = prep_parameters_arrays_CPM(all_mats,...
        all_covars, k, iterations);   
    else 
        [~, ~, int_pos_ntwrk, int_neg_ntwrk,...
        int_combined_ntwrk, slope_pos_ntwrk, slope_neg_ntwrk, ...
        slope_combined_ntwrk, slope_pos_covars, slope_neg_covars,...
        slope_combined_covars] = prep_parameters_arrays_CPM(all_mats,...
        all_covars, k, iterations);   
    end
    
    %% 4) Run CPM, evaluate model performance, extract selected edges + model parameters
    parfor i = 1:iterations
        fprintf('\n Running iteration # %6.3f\n',i);
         
        % Run CPM
         [behav_pred_pos_all(:,i), behav_pred_neg_all(:,i), ...
             behav_pred_combined_all(:,i), parameters_pos, parameters_neg,...
             parameters_combined, pos_mask, neg_mask, no_node, no_covars] = ...
             run_flexible_CPM(all_behav, all_mats, all_covars, k, ...
             thresh_type, thresh, adjust_stage);  %#ok<ASGLU>
        
        % Evaluate model performance
        [R_pos(i), P_pos(i), R_neg(i), P_neg(i), R_combined(i), P_combined(i),...
            rsq_pos(i), rsq_neg(i), rsq_combined(i), mae_pos(i), mae_neg(i),...
            mae_combined(i)] = evaluate_CPM(all_behav, behav_pred_pos_all(:,i),...
            behav_pred_neg_all(:,i), behav_pred_combined_all(:,i));
    
        % Extract model parameters
        if strcmp(adjust_stage, 'relate')
            [int_pos_ntwrk(i), int_neg_ntwrk(i), int_combined_ntwrk(i),...
                slope_pos_ntwrk(i), slope_neg_ntwrk(i), slope_combined_ntwrk(i),...
                ~, ~, ~] = extract_parameters_CPM(parameters_pos, parameters_neg,...
                parameters_combined, adjust_stage, no_covars); 
        else
            [int_pos_ntwrk(i), int_neg_ntwrk(i), int_combined_ntwrk(i),...
                slope_pos_ntwrk(i), slope_neg_ntwrk(i), slope_combined_ntwrk(i),...
                slope_pos_covars(i,:), slope_neg_covars(i,:), ...
                slope_combined_covars(i,:)] = ...
                extract_parameters_CPM(parameters_pos, parameters_neg,...
                parameters_combined,adjust_stage, no_covars); 
        end 
    
    % Store positive and negative edges masks - accounting for parallel loop
        c_pos_mask_all{i,1} = pos_mask;
        c_neg_mask_all{i,1} = neg_mask;  
    end
    
    % Add positive and negative edge masks
    pos_mask_all = [];
    neg_mask_all = [];
    for j = 1:iterations
        pos_mask_all = cat(3, c_pos_mask_all{j,1}, pos_mask_all);
        neg_mask_all = cat(3, c_neg_mask_all{j,1}, neg_mask_all);
    end
    
    %% 5) Extract selected edges 
    % set number of folds across all iterations of model
    k_all = k*iterations;
    
    % no_node not available after parfor loop so call again here
    no_node = size(all_mats,1);
    
    % get indices of edges in each network that are correlated with the
    % target variable in >= number of total folds specified by freq_thresh
    [pos_edges_all, neg_edges_all, pos_edges_thresh, neg_edges_thresh] = ...
        extract_edges_CPM(pos_mask_all, neg_mask_all, no_node, k_all, freq_thresh);
    
    %% 6) Store parameters, predicted values, edges, and results in structs and save data
    if strcmp(adjust_stage, 'relate')
        parameters = struct('int_pos_ntwrk', int_pos_ntwrk, 'int_neg_ntwrk', ...
        int_neg_ntwrk, 'int_combined_ntwrk', int_combined_ntwrk,...
        'slope_pos_ntwrk', slope_pos_ntwrk, 'slope_neg_ntwrk', slope_neg_ntwrk,...
        'slope_combined_ntwrk', slope_combined_ntwrk,...
        'pos_mask_all', pos_mask_all, 'neg_mask_all', neg_mask_all); 
    else
        parameters = struct('int_pos_ntwrk', int_pos_ntwrk, 'int_neg_ntwrk', ...
        int_neg_ntwrk, 'int_combined_ntwrk', int_combined_ntwrk,...
        'slope_pos_ntwrk', slope_pos_ntwrk, 'slope_neg_ntwrk', slope_neg_ntwrk,...
        'slope_combined_ntwrk', slope_combined_ntwrk, 'slope_pos_covars',...
        slope_pos_covars, 'slope_neg_covars', slope_neg_covars,...
        'slope_combined_covars', slope_combined_covars,...
        'pos_mask_all', pos_mask_all, 'neg_mask_all', neg_mask_all);
    end
    
    predictions = struct('pos_preds_all_folds', behav_pred_pos_all, ...
        'neg_preds_all_folds', behav_pred_neg_all, 'combined_preds_all_folds',...
        behav_pred_combined_all);
    
    edges = struct('pos_edges_all', pos_edges_all, 'neg_edges_all', neg_edges_all,...
        'pos_edges_thresh', pos_edges_thresh, 'neg_edges_thresh', neg_edges_thresh);
    
    results = struct('R_pos', R_pos, 'R_neg', R_neg, 'R_combined', R_combined,...
        'rsq_pos', rsq_pos, 'rsq_neg', rsq_neg, 'rsq_combined', rsq_combined,...
        'mae_pos', mae_pos, 'mae_neg', mae_neg', 'mae_combined', mae_combined);
    
    % Get mean parameters, predicted values, and results across iterations
    if iterations > 1
        mean_parameters = structfun(@mean, parameters, 'UniformOutput', false);
        mean_predictions = structfun(@(x) mean(x, 2), predictions, 'UniformOutput', false);
        mean_results = structfun(@mean, results, 'UniformOutput', false);
    end
    
    % save data
    variable_list = {'model_info', 'parameters', 'predictions', 'results', 'edges', ....
        'mean_parameters', 'mean_predictions', 'mean_results'};
    
    save_file = [model_info.output_path filesep model_info.model_name '.mat'];
    for variable_ix = 1:length(variable_list)
        if exist(variable_list{variable_ix})
            if exist(save_file)
                save(save_file, variable_list{variable_ix}, '-append');
            else
                save(save_file, variable_list{variable_ix});
            end
        end
    end
    
    %% 7) Create masks for visualisation
    % Check if restricted or computationally lesioned connectivity matrix used
    % and get indices of selected edges in full (i.e. 268 * 268 connectivity
    % matrix). 
    % NOTE: THIS CODE ASSUMES SHEN ATLAS PARCELLATION USED.
    if no_node < 268
                      
        % Read in csv file mapping node indices in restricted/lesioned
        % connectivity matrix to node indices in original connectivity matrix
        % 1st column is original index, 2nd column is new index
        node_indices = csvread('Y:\cogReserve\CPM_info\restricted_timeseries_node_indices.csv',1, 1);
    
        % map edge indices back to indices in full connectivity matrix
        [pos_edges_orig, neg_edges_orig] = get_original_edge_indices_CPM(...
            pos_edges_thresh, neg_edges_thresh, node_indices);
    
        % create binary edge masks (in 268*268 connectivity matrix)
        [pos_edge_mask, neg_edge_mask] = ...
        create_masks_CPM(pos_edges_orig, neg_edges_orig, 268);
    
    else  % if using a full timeseries 
        pos_edges_orig = pos_edges_thresh(:, [1 3 4]);
        neg_edges_orig = neg_edges_thresh(:, [1 3 4]);
        [pos_edge_mask, neg_edge_mask] = ...
        create_masks_CPM(pos_edges_orig, neg_edges_orig, 268);
    end
    
    % save binary edge masks for visualisation in bioimagesuite
    pos_mask_file = [model_info.output_path filesep model_info.model_name '_pos_mask.txt'];
    neg_mask_file = [model_info.output_path filesep model_info.model_name '_neg_mask.txt'];
    save(pos_mask_file, 'pos_edge_mask', '-ascii');
    save(neg_mask_file, 'neg_edge_mask', '-ascii');
    
    % save info on edges wrt to their position in the original 268*268
    % connectivity matrix
    orig_edge_ix = [model_info.output_path filesep model_info.model_name '_orig_edge_ix.mat'];
    save(orig_edge_ix, 'pos_edges_orig', 'pos_edge_mask', 'neg_edges_orig', 'neg_edge_mask');
    
    
    %% 8) Run permutation test and save
    if strcmp(run_permutation, 'yes')
        if iterations > 1
            [perm_p_pos, perm_p_neg, perm_p_combined] = ...
                CPM_permutation_test_parallelised(all_behav, all_mats, ...
                all_covars, k, thresh_type,thresh, adjust_stage,...
                mean_results.R_pos, mean_results.R_neg, mean_results.R_combined,...
                perm_iterations);
        else 
            [perm_p_pos, perm_p_neg, perm_p_combined] = ...
                CPM_permutation_test_parallelised(all_behav, all_mats, ...
                all_covars, k, thresh_type, thresh, adjust_stage, ...
                R_pos, R_neg, R_combined, perm_iterations);
        end
        
        % save permutation test results
        perm_results = struct('perm_p_pos', perm_p_pos, 'perm_p_neg', ...
            perm_p_neg, 'perm_p_combined', perm_p_combined); 
    
        save(save_file, 'perm_results', '-append');
    end

    All_pos(t,1) = mean_results.R_pos;
    All_neg(t,1) = mean_results.R_neg;
    All_combined(t,1) = mean_results.R_combined;
end
save('W:\CPM_FU3\ICV\output\DiffThreshold_both_Results.mat','All_pos',...
    'All_neg','All_combined','thresholds','k','iterations','adjust_stage');