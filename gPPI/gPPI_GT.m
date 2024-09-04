clear;
clc;

%% This script used to perform generalized PPI derived from Go trials

% Set up path
Para.rootdir = 'W:\All_IMAGEN_DATA\imaging\SST_spm_first_level_GT';  % SPM first level path
Para.opdir = 'W:\gPPI_Shen_GT';  % output path
Para.Maskdir = 'W:\Shen268_ROIs';% ROI mask path including 268 regions nfiti data
Para.subjfolder = dir([Para.rootdir '\00*']);
Para.ROInum = 268; % regions number
Para.Uu_Fail = [1,1,1]; % SPM.Sess.U(i).name{j} 'Fail_stop' [i,j,weight]
Para.Uu_Success = [2,1,1]; % SPM.Sess.U(i).name{j} 'Successful_stop'  [i,j,weight]
Para.Uu_Go = [3,1,1]; % SPM.Sess.U(i).name{j} 'Go trial'  [i,j,weight]
Data_miss.roi = [];
Data_miss.SPM = [];
Psy_Fail = [];


for  subji= 1 :length(subcode)
    tic
    %% Extract VOI 
    spm_file = fullfile(Para.rootdir,num2str(Para.subjfolder(subji),'%012d'),'SPM.mat');
    if exist(spm_file,'file') 
        load(spm_file);
        volumes = SPM.nscan;
        for roii = 1: Para.ROInum  % 268 ROIs
            % VOI job setting
            job.spmmat = {fullfile(Para.rootdir,num2str(Para.subjfolder(subji),'%012d'),'SPM.mat')};% SPM.mat path
            job.adjust = 1; % index of F-contrast
            job.session = 1;
            job.name = ['Shen268_' num2str(roii,'%.3d')];
            job.roi{1}.mask.image = {[Para.Maskdir filesep 'Shen_ROI_',num2str(roii,'%03d') ,'.nii,1']};% load mask.nii
            job.roi{1}.mask.threshold = 0.5;
            job.roi{2}.spm.spmmat = {''};
            job.roi{2}.spm.contrast = 1; % index of contrast
            job.roi{2}.spm.conjunction = 1;
            job.roi{2}.spm.threshdesc = 'none';
            job.roi{2}.spm.thresh = 1;
            job.roi{2}.spm.extent = 0;
            job.roi{2}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
            job.expression = 'i1&i2';
            try
                % Extract Y from VOI
                [Y, xY] = spm_regions_mask(job); 
                VOI(:,roii) = Y; 
                
               %% PPI-term calculation for each VOI
                % 'Fail_stop'PPI term
                Fail_PPI = spm_peb_ppi_Nosave(cell2mat(job.spmmat),'ppi',xY,Para.Uu_Fail,['Fail_gPPI_Shen_SST',num2str(roii)],0); % No PPI.mat saved
                % 'Success_stop'PPI term
                Success_PPI = spm_peb_ppi_Nosave(cell2mat(job.spmmat),'ppi',xY,Para.Uu_Success,['Success_gPPI_Shen_SST',num2str(roii)],0); % No PPI.mat saved
                % 'Go'PPI term
                Go_PPI = spm_peb_ppi_Nosave(cell2mat(job.spmmat),'ppi',xY,Para.Uu_Go,['Go_gPPI_Shen_SST',num2str(roii)],0); % No PPI.mat saved        
                
                PPI_Fail_all(:,roii) = Fail_PPI.ppi; % 'Fail_stop' ppi
                PPI_Success_all(:,roii) = Success_PPI.ppi;   % 'Success_stop' ppi  
                PPI_Go_all(:,roii) = Go_PPI.ppi;   % 'Success_stop' ppi  
                
                Psy_Fail = Fail_PPI.P;          % 'Fail_stop' psychological variables     
                Psy_Success = Success_PPI.P;    % 'Fail_Success' psychological variables 
                Psy_Go = Go_PPI.P;  
            catch 
                % If no voxels survise in mask & F-contrast.mask
                VOI(1:volumes,roii) = 0;  
                PPI_Fail_all(1:volumes,roii) = 0; % 'Fail_stop' ppi
                PPI_Success_all(1:volumes,roii) =0; % 'Success' ppi
                PPI_Go_all(1:volumes,roii) =0; % 'Success' ppi
                sprintf('Error: Subject %s fail to extract VOI %d!',num2str(Para.subjfolder(subji),'%012d'),roii)
                Data_miss.roi = [Data_miss.roi;strcat(num2str(subji,'%04d'), ',', num2str(Para.subjfolder(subji),'%012d'),',',num2str(roii,'%03d'))]; 
            end 
        end
                   
        if ~isempty(Psy_Fail)
        %% build PPI_GLM 
            for roii = 1: Para.ROInum
                if ~all(isnan(VOI(:,roii)))
                    PPI_matrix = [Psy_Fail Psy_Success Psy_Go VOI(:,roii) PPI_Fail_all(:,roii) PPI_Success_all(:,roii) PPI_Go_all(:,roii) ones(volumes,1)];
                    % (1)PsyVariable_Fail(2)PsyVariable_Success(3)VOI_BOLD(4)Fail_PPI_Term
                    % (5)Success_PPI_Term (6)constant
                    for roij = 1: Para.ROInum
                        if roij==roii
                            beta_PPI_SST(:,roii,roij) = zeros(8,1);
                            con_PPI_SST(roii,roij) = 0; % Success-Fail
                            con_PPI_Fail(roii,roij) = 0; % Fail stop
                            con_PPI_Success(roii,roij) = 0; % Success stop
                            con_PPI_Go(roii,roij) = 0; % Go trial
                        else
                            [b,bint] = regress(zscore(VOI(:,roij)),PPI_matrix);
                            beta_PPI_SST(:,roii,roij) = b;
                            con_PPI_SST(roii,roij) = [0 0 0 0 -1 1 0 0]*b; % Success-Fail
                            con_PPI_Fail(roii,roij) = [0 0 0 0 1 0 0 0]*b; % Fail stop
                            con_PPI_Success(roii,roij) = [0 0 0 0 0 1 0 0]*b; % Success stop
                            con_PPI_Go(roii,roij) = [0 0 0 0 0 0 1 0]*b; % Go trial
                        end
                    end  
                    
                % If no voxels survise in mask & F-contrast.mask
                else
                    for roij = 1: Para.ROInum
                        beta_PPI_SST(:,roii,roij) = zeros(8,1);
                        con_PPI_SST(roii,roij) = 0; % Success-Fail
                        con_PPI_Fail(roii,roij) = 0; % Fail
                        con_PPI_Success(roii,roij) = 0; % Success
                        con_PPI_Go(roii,roij) = 0; % Go

                    end  
                end
             end
           %% save file         
            mkdir(fullfile(Para.opdir,num2str(Para.subjfolder(subji),'%012d'))) ;
            save(fullfile(Para.opdir,num2str(Para.subjfolder(subji),'%012d'),'GLM_PPI_Shen_VOI268_SSTcon.mat'), 'beta_PPI_SST','con_PPI_SST','con_PPI_Fail','con_PPI_Success','con_PPI_Go');
            save(fullfile(Para.opdir,num2str(Para.subjfolder(subji),'%012d'),'PPI_Shen_VOI268.mat'), 'PPI_Fail_all','PPI_Success_all','VOI','Psy_Fail','Psy_Success','PPI_Go_all','Psy_Go');
            sprintf('Subject %s PPI_GLM have done!',num2str(Para.subjfolder(subji),'%012d'))   
        end  
    % No SPM.mat exist 
    else
         sprintf('Subject %s SPM missing!',num2str(Para.subjfolder(subji),'%012d'))
         Data_miss.SPM = [Data_miss.SPM;strcat(num2str(subji,'%04d'), ',', num2str(Para.subjfolder(subji),'%012d'))];
    end
    clearvars -except Para Data_miss
    Psy_Fail = [];

    toc
end 


