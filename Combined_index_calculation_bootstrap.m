%% Calculate the combined index. Produce visual tables from the grid results of RCs.
% Choose which dataset

addpath('C:\Users\user\Documents\Enrico\12_Consciousness_Connectomics\Bootstrap')
%% Adjust dimension of RC table, bootstrap 2
load bootstrap2_RCs.mat
%%
%RCs{33,13}=RCs{22,8};
%save 'bootstrap2_RCs.mat' RCs
%%
clc;
for numb=1:4
    series=[0, 1, 2, 8];
    bootstrap_id=series(numb);
    filename=sprintf('bootstrap%d_RCs.mat', bootstrap_id);
    disp(filename)
    load(filename)

% pca_var
xtick1=1:13;
xticklabels1=[0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95];
ytick1=1:33;

% Ratio of robust, weights_var, weights_var_real
ratio_of_robust=zeros(size(RCs,1),size(RCs,2));
weights_var=zeros(size(RCs,1),size(RCs,2));
weights_var_real=zeros(size(RCs,1),size(RCs,2));
weights_max_var_trait=zeros(size(RCs,1),size(RCs,2));

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

%  Combining the two measures into one using Quartiles
q=quantile(nonzeros(ratio_of_robust), [0.25 0.5 0.75 1]);
ratio_of_robust(q(3)<ratio_of_robust & ratio_of_robust<=q(4) & ratio_of_robust~=0)=[4];
ratio_of_robust(q(2)<ratio_of_robust & ratio_of_robust<=q(3) & ratio_of_robust~=0)=[3];
ratio_of_robust(q(1)<ratio_of_robust & ratio_of_robust<=q(2) & ratio_of_robust~=0)=[2];
ratio_of_robust(ratio_of_robust<=q(1) & ratio_of_robust~=0)=[1];

q=quantile(nonzeros(weights_var_real), [0.25 0.5 0.75 1]);
weights_var_real(q(3)<weights_var_real & weights_var_real<=q(4) & weights_var_real~=0)=[1];
weights_var_real(q(2)<weights_var_real & weights_var_real<=q(3) & weights_var_real~=0)=[2];
weights_var_real(q(1)<weights_var_real & weights_var_real<=q(2) & weights_var_real~=0)=[3];
weights_var_real(weights_var_real<=q(1) & weights_var_real~=0)=[4];

combined=ratio_of_robust+weights_var_real;

% Save files
path='C:\Users\user\Documents\Enrico\12_Consciousness_Connectomics\Bootstrap\Results';
fname = fullfile ( path, sprintf('combined_bootstrap%d.mat', bootstrap_id));
save (fname, 'combined')

fname = fullfile ( path, sprintf('ratio_of_rubust_bootstrap%d.mat', bootstrap_id));
save (fname, 'ratio_of_robust')

fname = fullfile ( path, sprintf('weights_var_bootstrap%d.mat', bootstrap_id));
save (fname, 'weights_var')

fname = fullfile ( path, sprintf('weights_var_real_bootstrap%d.mat', bootstrap_id));
save (fname, 'weights_var_real')

% Aggregate. Add up the bootstrap matrix results then divide to calculate mean

if numb==1
    ratio_of_robust_aggregate=ratio_of_robust;
    weights_var_aggregate=weights_var;
    weights_var_real_aggregate=weights_var_real;
end

ratio_of_robust_aggregate=ratio_of_robust_aggregate+ratio_of_robust;
weights_var_aggregate=weights_var_aggregate+weights_var;
weights_var_real_aggregate=weights_var_real_aggregate+weights_var_real;
end

ratio_of_robust_aggregate=ratio_of_robust_aggregate/numb;
weights_var_aggregate=weights_var_aggregate/numb;
weights_var_real_aggregate=weights_var_real_aggregate/numb;

combined_aggregate_bootstrap=ratio_of_robust_aggregate+weights_var_real_aggregate;

%% Visual plots. 
% Plot combined
figure; 
    imagesc(combined); colormap jet; colorbar; axis square;
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1);
    xlabel('Variance explained by PCA', 'fontsize', 10); ylabel('Number of IC', 'fontsize', 10);
    suptitle('Combined index, PCA/ICA grid, 8-most consistent model');
    %title('Sample dataset (Epilepsy + control)', 'fontsize', 12);
    title('Bootstrapped dataset 37/46 (sample)', 'fontsize', 12); 
    
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

%% Ratio of robust
close all; clc;
figure; 
    imagesc(ratio_of_robust); colormap hot; colorbar; axis square;
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1);
    xlabel('Variance explained by PCA', 'fontsize', 10); ylabel('Number of IC', 'fontsize', 10);
    suptitle('Ratio of robust components obtained from each model');
    %title('Sample dataset (healthy control + epilepsy)', 'fontsize', 12);   
    title('Bootstrapped dataset 37/46 (sample)', 'fontsize', 12); 
%% Variance of weights for each model (log)
figure; 
    imagesc((weights_var)); colormap jet; colorbar; axis square;
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1);
    xlabel('Variance explained by PCA', 'fontsize', 10); ylabel('Number of IC', 'fontsize', 10);
    suptitle('Variance of weights for each model')
    title('Bootstrapped dataset 37/46 (sample)', 'fontsize', 12); 

%% Relative variance of weights for each model 
close all;
figure; 
    imagesc((weights_var_real), [0,0.3]); myColorMap = flipud(hot); myColorMap(1,:) = 1; colormap(myColorMap);
    colorbar; axis square;
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1);
    xlabel('Variance explained by PCA', 'fontsize', 10); ylabel('Number of IC', 'fontsize', 10);
    suptitle('Relative variance of weights for each model  ');
    %title('Sample dataset (healthy control + epilepsy)', 'fontsize', 12);  
    title('Bootstrapped dataset 37/46 (sample)', 'fontsize', 12); 

