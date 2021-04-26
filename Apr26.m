%% Mihály Gayer, University of Geneva, 2021.02.02 Master thesis project, 
% under the supervision of Dr Enrico Amico and Jonathan Wirschi
% 1. load connICA matrix (obtained from connICA_and_FC_construction.mat file)
% 2. Configure ICA setup
% 3. To obtain the RC: run_global_robust_connICA (function is saved separately) 
% 4. Plot the robust components and the weights mean, median with statistical tests 
close all; clearvars; clc; 
%% 1. load connICA matrix
addpath('C:\Users\user\Documents\Enrico\12_Consciousness_Connectomics')
%% - Jonathan's dataset
configs.which_dataset='control'; % 1=Jonathan   2 = class
load connICA_matrix_control.mat 
load FC_control.mat

%% Class - Dataset
configs.which_dataset='class'; % 1=Jonathan   2 = class
load connICA_matrix_class.mat
load FC_class.mat

%% Compiled 2 healthy datasets
configs.which_dataset='2healthy'; % 1=Jonathan   2 = class
load connICA_matrix_2healthy.mat
load FC.mat

%% Epilepsy 
% rTLE patient 6 could be outlier
configs.which_dataset='epilepsy'; % 
load connICA_matrix_epilepsy.mat
%load FC_epilepsy.mat
%% Sample = Epilepsy + Jonathan control
configs.which_dataset='sample';
load connICA_matrix_sample

%% Sample = Epilepsy + Jonathan control
configs.which_dataset='class_epilepsy';
load connICA_matrix_class_epilepsy.mat

%% All = Epilepsy + Control + Class
configs.which_dataset='all';
load connICA_matrix_all

%% 2. ICA Set up and configuration
load('eeg-fmri_connectomes_destrieux_scrubbed_256Ch3T.mat')
load('atlas_labels.mat')
load('aparc_a2009_yeoRS7_148reg_eeg_nosubc.mat') % from Jonathan's email, 2021.02.01
%set(0,'DefaultFigureWindowStyle','docked'); % to dock figures within command window
set(0,'DefaultFigureWindowStyle','normal'); % to go back to normal mode
% ConnICA params
configs.numRegions = 148; % aparc2009 parcellation
configs.mask = triu(true(configs.numRegions,configs.numRegions),1); % upper triangular are ones in 148x148 matrix
configs.numConn = size(connICA_matrix,1);%number of connectomes (subjects or sessions)
configs.numEdges = nnz(configs.mask);
configs.epsilon = 0.0001; configs.maxFinetune=100; configs.bootstrap=1;
addpath(genpath('FastICA_25'));  
%add FastICA folder to path. Please check folder_name and update if needed.
% ICA parameters
configs.maxNumIterations=300;
configs.numRuns = 100;           % run FastICA several time for robust comps extraction  
configs.numRuns_max=300;          % maximum number of times fastICA runs per 1 set-up. Aim is to obtain numRuns amount. numRuns should be smaller than numRuns_max ALWAYS.

configs.minFreq = 0.75;          % minimum frequency required 
configs.corrMin = 0.75;          % correlation between FC_traits
configs.corrMin_A = 0.75;        % correlation between subject weights, set it to 0 if bootstrapping subjects across numRuns 

%% 3. run_global_robust_connICA
%% INPUT CONNICA MATRIX
connICA_matrix=connICA_matrix_sample;
%% Double check PCA variance structure and apply the corresponding last.Eig
configs.lastEig=35;
configs.numOfIC=30;
%%
clc; clearvars RC RC_Index
[RC, RC_Index]=run_global_robust_connICA(connICA_matrix, configs);
clc;
disp('done')
disp(length(RC))

%%
%save 'RC_sample_31_13.mat' RC
%%
load RC_sample_28_13.mat
%%
load RC_sample_31_13.mat


%% 4. Plot the robust components and the weights mean, median with statistical tests 
%% Plot one RC 
t=4
figure; 
    imagesc(RC(t).matrix(yeoOrder,yeoOrder),[-3,3]); colormap jet; colorbar; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); xlabel('regions'); ylabel('regions');
    title(sprintf('connICA comp %d',t));
    
%% Plot All RC traits
hold off
tiledlayout(2,ceil(length(RC)/2))
for t = 1:length(RC)-1
  nexttile
    imagesc(RC(t).matrix(yeoOrder,yeoOrder),[-3,3]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title(sprintf('IC %d',t)); hold on
end
nexttile
    imagesc(RC(length(RC)).matrix(yeoOrder,yeoOrder),[-3,3]); colormap jet; axis square; colorbar
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title(sprintf('IC %d',length(RC))); hold on
%% Plot mean weights of RC
hold off
plt=tiledlayout(2,ceil(length(RC)/2));
for t = 1:length(RC)
 nexttile
    plot(1:length(RC(t).weights_mean), RC(t).weights_mean);
    xline(21, 'red'); 
    title(sprintf('IC %d', t));
end
ylabel(plt, 'mean weights');
xlabel(plt, 'subject id, healthy control: 1-21;  epileptic patients: 22-46');
hold off

%% Median - t test
ps_ttest=zeros(length(RC),1);
hold off
plt=(tiledlayout(2,ceil(length(RC)/2)))
for t = 1:length(RC)
    mdn=nanmedian(RC(t).weights,2);
    [h,p_ttest]=ttest2(nanmedian(RC(t).weights(1:21,:),2),nanmedian(RC(t).weights(22:46,:),2), 'Vartype','unequal');
    ps_ttest(t)=round(p_ttest,2);
    
 nexttile
    plot(1:46, mdn);
    xline(21, 'red'); 
    title(['IC ',num2str(t),'; t-test p value: ',num2str(ps_ttest(t))]);
end
ylabel(plt, 'meadian weights');
xlabel(plt, 'subject id,  epileptic patients: 22-46');
hold off
%% Median - Wilcoxon test
ps_wilcoxon=zeros(length(RC),1);
hold off
tiledlayout(2,ceil(length(RC)/2))
for t = 1:length(RC)
    mdn=nanmedian(RC(t).weights,2);
  
    [p_wilcoxon,h]=ranksum(nanmedian(RC(t).weights(1:21,:),2),nanmedian(RC(t).weights(22:46,:),2));
    ps_wilcoxon(t)=round(p_wilcoxon,2);
    
 nexttile
    plot(1:46, mdn); ylabel('median weights'); xlabel('subject id');
    xline(21, 'red'); 
    title(['IC ',num2str(t),'; Wilcoxon p value: ',num2str(ps_wilcoxon(t))]);
end
hold off
%% Boxplot
boxplot(RC(7).weights'); xlabel('subject id, 1-21 healthy subjects'); ylabel('weights');

%%
clc;
t=5;
[h,p]=ttest2(nanmedian(RC(t).weights(1:21,:),2),nanmedian(RC(t).weights(22:46,:),2), 'Vartype','unequal')
%%
ttest2(RC(t).weights(1:21),RC(t).weights(22:46))