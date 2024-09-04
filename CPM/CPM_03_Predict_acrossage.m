clear;
clc;

% This codes used models identified from previous timepoint to predict
% the following timepoint
train_group = 'BL';  % BL (age 14), FU2 (age 19), FU3 (age 23)
test_group = ' FU2'; % FU2 (age 19), FU3 (age 23)

% covariable "pd (Puberty development scale)" index 4 for BL group(age 14);
% if train_group is FU2 or FU3, id_pd = [];
id_pd = 4;  

% output file
output_path = fullfile('W:\Results_Shen_GT\General_results\ICV');
mkdir(output_path);
output_file = fullfile(output_path,[train_group,'_predict_',test_group,'.mat']);

% train mats,covars,behav
load(fullfile('W:\Results_Shen_GT\ICV',train_group,'cpm_mat\cpm_predictor_variables_all_ppts.mat'),'cpm_predictors');
train_mats = cpm_predictors;
load(fullfile('W:\Results_Shen_GT\ICV',train_group,'cpm_mat\cpm_target_variable.mat'), 'cpm_target');
train_behav = cpm_target;
train_covars = readtable(fullfile('W:\Yihe\Both_Results_Shen_GT\ICV',train_group,'cpm_mat\site_covariates.csv'));
train_covars = train_covars{:,2:end};

% testing mats,covars,behav
load(fullfile('W:\Results_Shen_GT\ICV',test_group,'cpm_mat\cpm_predictor_variables_all_ppts.mat'),'cpm_predictors');
test_mats = cpm_predictors;
load(fullfile('W:\Results_Shen_GT\ICV',test_group,'cpm_mat\cpm_target_variable.mat'), 'cpm_target');
test_behav = cpm_target; 
test_covars = readtable(fullfile('W:\Yihe\Results_Shen_GT\ICV',test_group,'cpm_mat\site_covariates.csv'));
test_covars.pds = []; 
test_covars = test_covars{:,2:end};
n = length(test_behav);

% CPM parameters
adjust_stage = 'both';
thresh = 0.01;
thresh_type = 'p-value';

% apply the models from previous timepoint (train data) to the following
% timepoint (test data)
[R_pos,R_neg,R_combined,pred_pos,pred_neg,pred_combined,pos_mask,neg_mask]...
    = General_cpm(train_group,test_group,train_mats,train_behav,train_covars,...
    test_mats,test_behav,test_covars,adjust_stage,thresh_type,thresh,id_pd);

model_info = struct('train_group',train_group,'test_group',test_group, ...
    'adjust_stage', adjust_stage, 'thresh', thresh,'thresh_type',thresh_type);


%% save binary edge masks for visualisation in bioimagesuite
pos_mask_file = [output_path filesep train_group,'_predict_',test_group, '_pos_mask.txt'];
neg_mask_file = [output_path filesep train_group,'_predict_',test_group, '_neg_mask.txt'];
save(pos_mask_file, 'pos_mask', '-ascii');
save(neg_mask_file, 'neg_mask', '-ascii');

%% permutation test
perm_iterations = 1000;
random_pos_r = zeros(perm_iterations,1);
random_neg_r = zeros(perm_iterations,1);
random_combined_r = zeros(perm_iterations,1);

for i = 1:perm_iterations
    random_pos_behav = pred_pos(randperm(n));
    random_neg_behav = pred_neg(randperm(n));
    random_combined_behav = pred_combined(randperm(n));

    random_pos_r(i) = corr(test_behav,random_pos_behav);
    random_neg_r(i) = corr(test_behav,random_neg_behav);
    random_combined_r(i) = corr(test_behav,random_combined_behav);

end

% Calculate permuted p-value for positive network strength model
perm_p_pos = (sum(random_pos_r>=R_pos))/perm_iterations;

% Calculate permuted p-value for negative network strength model 
perm_p_neg = (sum(random_neg_r>=R_neg))/perm_iterations;

% Calculate permuted p-value for combined network strength model 
perm_p_combined = (sum(random_combined_r>=R_combined))/perm_iterations;  

save(output_file,'R_pos','R_neg','R_combined','model_info','perm_p_pos','perm_p_neg','perm_p_combined',...
    'pred_pos','pred_neg','pred_combined','test_behav');
