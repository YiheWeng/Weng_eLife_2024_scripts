function [R_pos,R_neg,R_combined,pred_pos,pred_neg,pred_combined,pos_mask,neg_mask]...
    =General_cpm(train_group,test_group,train_mats,train_behav,train_covars,...
    test_mats,test_behav,test_covars,adjust_stage,thresh_type,thresh,id_pd)
    %% CPM analysis
    no_node = size(train_mats,1);
    test_no_covars = size(test_covars,2);
    train_no_covars = size(train_covars,2);

    train_vcts = reshape(train_mats,[],size(train_mats,3));

    if strcmp(adjust_stage, 'relate') | strcmp(adjust_stage, 'both')
        [r_mat, p_mat] = CPM_fs_relate_partial(train_vcts, train_behav, ...
            train_covars, no_node);
    else
        [r_mat, p_mat] = CPM_fs_relate(train_vcts, train_behav, no_node);
    end

    % feature selection - select edges (Step 4 - Shen et al. 2017)
    if strcmp(thresh_type, 'p-value')
        [pos_mask, neg_mask] = CPM_fs_select_pvalue(r_mat, p_mat,...
            thresh, no_node);

    else strcmp(thresh_type, 'sparsity');
        [pos_mask, neg_mask] = fs_select_sparsity_CPM(r_mat, p_mat,...
            thresh, no_node);  
    end

    %% network strength
    for ss = 1:size(train_mats,3)
        train_sumpos(ss) = sum(sum(train_mats(:,:,ss).*pos_mask))/2;
        train_sumneg(ss) = sum(sum(train_mats(:,:,ss).*neg_mask))/2;
    end

    % Calculate combined network strength
    train_sumcombined = train_sumpos - train_sumneg;

    % fit model on training set (Step 6 - Shen et al. 2017)
    % include covars in multiple regression if specified - otherwise use 
    % simple linear regression
    if strcmp(adjust_stage, 'fit') | strcmp(adjust_stage, 'both')
        [fit_pos, fit_neg, fit_combined] = ...
        CPM_fit_model(train_behav, train_covars, train_no_covars, ...
        train_sumpos', train_sumneg', train_sumcombined'); 
        % delete BL pds slope to predict FU2/FU3
        if strcmp(train_group, 'BL')
            fit_pos(id_pd+2) = [];
            fit_neg(id_pd+2) = [];
            fit_combined(id_pd+2) = [];
        end
        if strcmp(test_group, 'Stratify')
            fit_pos(6:end) = [];
            fit_neg(6:end) = [];
            fit_combined(6:end) = [];
        end         
    else
        [fit_pos, fit_neg, fit_combined] = ...
        CPM_fit_model(train_behav, [], 0, ...
        train_sumpos', train_sumneg', train_sumcombined');
    end
    
    % apply model to test set (Step 7 - Shen et al., 2017)
    % account for covars in model application if adjusted for in Step 6
    if strcmp(adjust_stage, 'fit') | strcmp(adjust_stage, 'both')
        [pred_pos, pred_neg, pred_combined] = ...
            CPM_apply_model(test_mats, test_covars, test_no_covars, ...
            pos_mask, neg_mask, fit_pos, fit_neg, fit_combined);
    else
        [pred_pos, pred_neg, pred_combined] = ...
            CPM_apply_model(test_mats, [], 0, ...
            pos_mask, neg_mask, fit_pos, fit_neg, fit_combined);
    end

    % Pearson's correlations
    [R_pos, P_pos] = corr(pred_pos,test_behav);
    [R_neg, P_neg] = corr(pred_neg,test_behav);
    [R_combined, P_combined] = corr(pred_combined,test_behav);
end

