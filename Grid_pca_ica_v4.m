%% After setting up the environment for the ICA do a grid ICA PCA result simulation.
%clc; clearvars;
%% Dataset: epilepsy + healthy control (Zenodo)

configs.which_dataset='sample'; % 
load connICA_matrix_sample.mat
connICA_matrix=connICA_matrix_sample;

%% PCA goes from 2:31 for 36%-99%
%% Loop goes for long

%% Repeat the algorithm for those cases where there were not enough runs to converge

configs.numRuns_max=1000; 
for i=1:40
    for k=1:i        
        if isequal(RCs_sample{k,i}, 'numRuns_max')
            configs.lastEig=i;
            configs.numOfIC=k;
            RCs_sample{k,i}=run_global_robust_connICA(connICA_matrix, configs); 
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bootstrap
clc;
addpath('C:\Users\user\Documents\Enrico\12_Consciousness_Connectomics')
load connICA_matrix_control.mat
load connICA_matrix_epilepsy.mat
% Take 17 from the control group (21), take 20 from the epileptic group (25)
%connICA_matrix=[connICA_matrix_control(randperm(21,17),:); connICA_matrix_epilepsy(randperm(25,20),:)];

configs.maxNumIterations=300;
configs.numRuns_max=800;
configs.numRuns = 100;   
%% 21 control, 25 epileptic patients 21:25 0.84 ~= 17:20

for j=1:1 % specify number of times to obtain bootstrapped RCs table
% connICA_matrix=[connICA_matrix_control(randperm(21,17),:); connICA_matrix_epilepsy(randperm(25,20),:)];
% name=sprintf('bootstrap%d.mat',j);
% configs.which_dataset=name; 
% configs.connICA_matrix=connICA_matrix;
% 
% % Each bootstrapped connICA has a different PCA structure, so recalculate
% %   .. how many eiganvectors we want to keep for the desired % of variance explained
% numPCAComps = size(connICA_matrix,1);
%     [~, ~, latent] = pca(connICA_matrix','NumComponents',numPCAComps); 
%     variance = cumsum(latent)./sum(latent); 
%     variance = variance(1:numPCAComps);
% pca_var=[0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95];
% pca_eigenvectors=zeros(1,numel(pca_var));
% 
% for i_pca=1:numel(pca_var)
%     pca_eigenvectors(1,i_pca)=find(variance>pca_var(i_pca),1);
% end
% 
% RCs=cell(pca_eigenvectors(numel(pca_eigenvectors)),numel(pca_var));

for i=11:13%numel(pca_var)
    lastEig=pca_eigenvectors(i);
    for k=2:lastEig
        configs.lastEig=lastEig;
        configs.numOfIC=k;
        [RC, RC_Index]=run_global_robust_connICA(connICA_matrix, configs);
        RCs{k,i}=RC;
    end
        % save intermediary results, after finishing a column (pca)
        filename=sprintf('bootstrap%d_RCs_pca_%d.mat',j, i ); 
        save(['C:\Users\user\Documents\Enrico\12_Consciousness_Connectomics\Bootstrap\Intermediary results\' filename],'RCs');

end 
filename=sprintf('bootstrap%d_RCs.mat',j);
save(['C:\Users\user\Documents\Enrico\12_Consciousness_Connectomics\Bootstrap\' filename],'RCs');

end
    

clc; disp('done')
%%

%% Produce visual tables from the grid results of RCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Choose which dataset
load RCs_sample_40.mat
%%
load pca_result_sample.mat
%%
xticklabels1=pca_result(1,:);
xtick1=1:40;
ytick1=1:40;
%%
load RCs2_v1.mat
RCs=RCs2;
load pca_result_class.mat
xticklabels2=pca_result_class(1,:);
xtick2=1:20;
ytick2=1:20;
%%
load RCs3_v1_26.mat
RCs=RCs3;
RCs=RCs(1:26,1:26);
load pca_result_2healthy.mat % double check the file name
xticklabels3=pca_result(1,1:26);
xtick3=1:26;
ytick3=1:26;
%% Ratio of robust, weights_var, weights_var_real
clc;
%RCs=RCs_sample;

ratio_of_robust=zeros(size(RCs,1));
weights_var=zeros(size(RCs,1));
weights_var_real=zeros(size(RCs,1));
weights_max_var_trait=zeros(size(RCs,1));
for i=1:size(RCs,1)
    for k=1:size(RCs,2)
        if isstruct(RCs{i,k})
            ratio_of_robust(i,k)=round( size(RCs{i,k},2) / RCs{i,k}(1).configs.numOfIC,2);
            weights_var(i,k)=mean([RCs{i,k}.weights_var]);
            weights_var_real(i,k)=mean([RCs{i,k}.weights_var_real]);
            weights_max_var_trait(i,k)=min([RCs{i,k}.inf_trait]);
        end
    end
end
%% Observing the relative variance weights matrix
% 97% of the relative weights are in the range of [0, 0.05], in general the
% weights differ by less than 5% from the mean
disp(length(weights_var_real(weights_var_real~=0)))% 690 non-zero obs 
disp(max(weights_var_real(:)))                     % 30.36 max value
disp(length(weights_var_real(weights_var_real>5))) % 2 obs larger than 5
disp(length(weights_var_real(weights_var_real>3))) % 3 obs larger than 3
disp(length(weights_var_real(weights_var_real>2))) % 5 obs larger than 2
disp(length(weights_var_real(weights_var_real>1))) % 8 obs larger than 1 (1-99%)
disp(length(weights_var_real(weights_var_real>0.5)))% 25 obs larger than 0.5
%% 
close all; clc
boxplot(weights_var_real(weights_var_real~=0 & weights_var_real<1))
suptitle('Relative variance of weights for each model')
title('Sample dataset (healthy control + epilepsy)', 'fontsize', 12);   
%%


%% Tick labels: load pca_result_sample.mat object, obtained from PCA_connICA.m script separately
ytick=1:size(RCs,1);
% % Sample
% xtick=[1, 6, 9, 12, 15, 18, 22, 25, 30, 35, 40];
% xticklabels1={'35%', '50%', '55%', '60%', '65%', '70%', '75%', '80%', '85%', '90%', '95%'};

% Bootstrap37_1
xtick=[1, 5, 7, 10, 12, 15, 18, 21, 24, 28, 32];
xticklabels1={'35%', '50%', '55%', '60%', '65%', '70%', '75%', '80%', '85%', '90%', '95%'};

%% Ratio of robust
close all; clc;
figure; 
    imagesc(ratio_of_robust); colormap jet; colorbar; axis square;
    set(gca,'xtick',xtick, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick);
    xlabel('Variance explained by PCA', 'fontsize', 10); ylabel('Number of IC', 'fontsize', 10);
    suptitle('Ratio of robust components obtained from each model');
    %title('Sample dataset (healthy control + epilepsy)', 'fontsize', 12);   
    title('Bootstrapped dataset 37/46 (sample)', 'fontsize', 12); 
%% Variance of weights for each model (log)
figure; 
    imagesc((weights_var)); colormap jet; colorbar; axis square;
    set(gca,'xtick',xtick, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick);
    xlabel('Variance explained by PCA', 'fontsize', 10); ylabel('Number of IC', 'fontsize', 10);
    suptitle('Variance of weights for each model')
    title('Bootstrapped dataset 37/46 (sample)', 'fontsize', 12); 
%% Relative variance of weights for each model (log)
figure; 
    imagesc(log(weights_var_real)); colormap jet; colorbar; axis square;
    set(gca,'xtick',xtick, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick);
    xlabel('Variance explained by PCA', 'fontsize', 10); ylabel('Number of IC', 'fontsize', 10);
    suptitle('Relative variance of weights for each model  ');
    %title('Sample dataset (healthy control + epilepsy)', 'fontsize', 12);  
    title('Bootstrapped dataset 37/46 (sample)', 'fontsize', 12); 
%% SAVE
% save 'ratio_of_robust_sample.mat' ratio_of_robust
% save 'weights_var_real_sample.mat' weights_var_real
load ratio_of_robust_sample.mat
load weights_var_real_sample.mat

% save 'ratio_of_robust_bootstrapped37_1.mat' ratio_of_robust
% save 'weights_var_real_bootstrapped37_1.mat' weights_var_real
% save 'weights_var_bootstrapped37_1.mat' weights_var


%%%%%%%%%%%%%%%%%%%%%
%%  04.07 Combining the two measures into one

%% PREVIOUS DATA
clc;
ratio_of_robust1=load('ratio_of_robust1.mat');
ratio_of_robust1=ratio_of_robust1.ratio_of_robust;

ratio_of_robust2=load('ratio_of_robust2.mat');
ratio_of_robust2=ratio_of_robust2.ratio_of_robust;

ratio_of_robust3=load('ratio_of_robust3.mat');
ratio_of_robust3=ratio_of_robust3.ratio_of_robust;
%
weights_var_real1=load('weights_var_real1.mat');
weights_var_real1=weights_var_real1.weights_var_real;

weights_var_real2=load('weights_var_real2.mat');
weights_var_real2=weights_var_real2.weights_var_real;

weights_var_real3=load('weights_var_real3.mat');
weights_var_real3=weights_var_real3.weights_var_real;
%% PREVIOUS DATA
entries=nonzeros(weights_var_real1); 
scatter(1:97,entries([1:90 92:96 98:99])) % 91 and 97 are outliers

q=quantile(nonzeros(ratio_of_robust1), [0.25 0.5 0.75 1]);
ratio_of_robust1(q(3)<ratio_of_robust1 & ratio_of_robust1<=q(4) & ratio_of_robust1~=0)=[4];
ratio_of_robust1(q(2)<ratio_of_robust1 & ratio_of_robust1<=q(3) & ratio_of_robust1~=0)=[3];
ratio_of_robust1(q(1)<ratio_of_robust1 & ratio_of_robust1<=q(2) & ratio_of_robust1~=0)=[2];
ratio_of_robust1(ratio_of_robust1<=q(1) & ratio_of_robust1~=0)=[1];

q=quantile(nonzeros(ratio_of_robust2), [0.25 0.5 0.75 1]);
ratio_of_robust2(q(3)<ratio_of_robust2 & ratio_of_robust2<=q(4) & ratio_of_robust2~=0)=[4];
ratio_of_robust2(q(2)<ratio_of_robust2 & ratio_of_robust2<=q(3) & ratio_of_robust2~=0)=[3];
ratio_of_robust2(q(1)<ratio_of_robust2 & ratio_of_robust2<=q(2) & ratio_of_robust2~=0)=[2];
ratio_of_robust2(ratio_of_robust2<=q(1) & ratio_of_robust2~=0)=[1];

q=quantile(nonzeros(ratio_of_robust3), [0.25 0.5 0.75 1]);
ratio_of_robust3(q(3)<ratio_of_robust3 & ratio_of_robust3<=q(4) & ratio_of_robust3~=0)=[4];
ratio_of_robust3(q(2)<ratio_of_robust3 & ratio_of_robust3<=q(3) & ratio_of_robust3~=0)=[3];
ratio_of_robust3(q(1)<ratio_of_robust3 & ratio_of_robust3<=q(2) & ratio_of_robust3~=0)=[2];
ratio_of_robust3(ratio_of_robust3<=q(1) & ratio_of_robust3~=0)=[1];

q=quantile(nonzeros(weights_var_real1), [0.25 0.5 0.75 1]);
weights_var_real1(q(3)<weights_var_real1 & weights_var_real1<=q(4) & weights_var_real1~=0)=[1];
weights_var_real1(q(2)<weights_var_real1 & weights_var_real1<=q(3) & weights_var_real1~=0)=[2];
weights_var_real1(q(1)<weights_var_real1 & weights_var_real1<=q(2) & weights_var_real1~=0)=[3];
weights_var_real1(weights_var_real1<=q(1) & weights_var_real1~=0)=[4];

q=quantile(nonzeros(weights_var_real2), [0.25 0.5 0.75 1]);
weights_var_real2(q(3)<weights_var_real2 & weights_var_real2<=q(4) & weights_var_real2~=0)=[1];
weights_var_real2(q(2)<weights_var_real2 & weights_var_real2<=q(3) & weights_var_real2~=0)=[2];
weights_var_real2(q(1)<weights_var_real2 & weights_var_real2<=q(2) & weights_var_real2~=0)=[3];
weights_var_real2(weights_var_real2<=q(1) & weights_var_real2~=0)=[4];

q=quantile(nonzeros(weights_var_real3), [0.25 0.5 0.75 1]);
weights_var_real3(q(3)<weights_var_real3 & weights_var_real3<=q(4) & weights_var_real3~=0)=[1];
weights_var_real3(q(2)<weights_var_real3 & weights_var_real3<=q(3) & weights_var_real3~=0)=[2];
weights_var_real3(q(1)<weights_var_real3 & weights_var_real3<=q(2) & weights_var_real3~=0)=[3];
weights_var_real3(weights_var_real3<=q(1) & weights_var_real3~=0)=[4];

combined1=ratio_of_robust1+weights_var_real1
combined2=ratio_of_robust2+weights_var_real2
combined3=ratio_of_robust3+weights_var_real3

%save 'combined1.mat' combined1
%save 'combined2.mat' combined2
%save 'combined3.mat' combined3

%% Sample data creat new ratio_of_robust1 and weights_var_real1 matrices
ratio_of_robust1=ratio_of_robust;
weights_var_real1=weights_var_real;

q=quantile(nonzeros(ratio_of_robust1), [0.25 0.5 0.75 1]);
ratio_of_robust1(q(3)<ratio_of_robust1 & ratio_of_robust1<=q(4) & ratio_of_robust1~=0)=[4];
ratio_of_robust1(q(2)<ratio_of_robust1 & ratio_of_robust1<=q(3) & ratio_of_robust1~=0)=[3];
ratio_of_robust1(q(1)<ratio_of_robust1 & ratio_of_robust1<=q(2) & ratio_of_robust1~=0)=[2];
ratio_of_robust1(ratio_of_robust1<=q(1) & ratio_of_robust1~=0)=[1];

q=quantile(nonzeros(weights_var_real1), [0.25 0.5 0.75 1]);
weights_var_real1(q(3)<weights_var_real1 & weights_var_real1<=q(4) & weights_var_real1~=0)=[1];
weights_var_real1(q(2)<weights_var_real1 & weights_var_real1<=q(3) & weights_var_real1~=0)=[2];
weights_var_real1(q(1)<weights_var_real1 & weights_var_real1<=q(2) & weights_var_real1~=0)=[3];
weights_var_real1(weights_var_real1<=q(1) & weights_var_real1~=0)=[4];

combined=ratio_of_robust1+weights_var_real1;
%% Save
save 'combined_sample.mat' combined
%save 'combined_bootstrapped37_1.mat' combined
%% Entire RCs_sample
figure; 
    imagesc(combined); colormap jet; colorbar; axis square;
    set(gca,'xtick',xtick, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick);
    xlabel('Variance explained by PCA', 'fontsize', 10); ylabel('Number of IC', 'fontsize', 10);
    suptitle('Combined index, PCA/ICA grid, 8-most consistent model');
    %title('Sample dataset (Epilepsy + control)', 'fontsize', 12);
    title('Bootstrapped dataset 37/46 (sample)', 'fontsize', 12); 

%% 75% + PCA 22:40 for sample, 18:37 for bootstrapped37_1
close all;
figure; 
    imagesc(combined(:, 22:40)); colormap jet; colorbar; axis square;
    set(gca,'xtick',[1, 3, 8, 13, 18], 'xticklabel', xticklabels1(7:11), 'fontsize',6); set(gca,'ytick',ytick1);
    xlabel('Variance explained by PCA', 'fontsize', 10); ylabel('Number of IC', 'fontsize', 10);
    suptitle('Combined index, PCA/ICA grid, 8-most consistent model');
    title('Sample dataset (Epilepsy + control)', 'fontsize', 12);

    
%% 75% + PCA 22:40 for sample, 18:37 for bootstrapped37_1
close all;
figure; 
    imagesc(combined(:, 18:32)); colormap jet; colorbar; axis square;
    set(gca,'xtick',[1, 4, 7, 11, 15], 'xticklabel', xticklabels1(7:11), 'fontsize',6); %set(gca,'ytick',ytick1);
    xlabel('Variance explained by PCA', 'fontsize', 10); ylabel('Number of IC', 'fontsize', 10);
    suptitle('Combined index, PCA/ICA grid, 8-most consistent model');
    %title('Sample dataset (Epilepsy + control)', 'fontsize', 12);
    title('Bootstrapped dataset 37/46 (sample)', 'fontsize', 12); 
    
    %% Previous data
figure; 
    imagesc(combined2); colormap jet; colorbar; axis square;
    set(gca,'xtick',xtick2, 'xticklabel', xticklabels2, 'fontsize',6); set(gca,'ytick',ytick2);
    xlabel('Variance explained by PCA', 'fontsize', 10); ylabel('Number of IC', 'fontsize', 10);
    title('Combined index, PCA/ICA grid Dataset 2 (Class) ', 'fontsize', 12)

%%
figure; 
    imagesc(combined3); colormap jet; colorbar; axis square;
    set(gca,'xtick',xtick3, 'xticklabel', xticklabels3, 'fontsize',6); set(gca,'ytick',ytick3);
    xlabel('Variance explained by PCA', 'fontsize', 10); ylabel('Number of IC', 'fontsize', 10);
    title('Combined index, PCA/ICA grid Dataset 3 (Healthy control all.) ', 'fontsize', 12)
