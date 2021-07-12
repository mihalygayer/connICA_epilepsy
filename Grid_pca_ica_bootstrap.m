%% After setting up the environment for the ICA do a grid ICA PCA result simulation.
% Bootstrap without replacement

% 1. load connICA matrix for epilepsy and control, create an example sample with replacement
% 2. setup configs for PCA and ICA
% 3. Run a PCA experimental calculation on 100 samples, how the the number of dimensions (last eigenvector)
%    required for a 95% PCA variance reduction. -> Grid dimension changes with replacement

% 4. Run the grid with replacement resampling. One RCs grid is with one sample
% 4.a Parfor at the number of Independent components level (rows in the grid)
%     Takes 440 minutes for one grid (one resample)

% 4.b Parfor at both number of IC and PCA - create a linear grid vector(s) and use on parfor loop only
%     Takes 30% less time than the other two options.

% 4.c Regular for loop took 480 minutes for one grid.

%% 1. Load connICA matrices, set up one example 
clc; 
clearvars;
addpath('C:\Users\user\Documents\Enrico\12_Consciousness_Connectomics')

load connICA_matrix_control.mat
load connICA_matrix_epilepsy.mat
% 21 control, 25 epileptic patients in the original sample
sample_vector_control=datasample(1:21,21,'Replace',true);
sample_vector_epilepsy=datasample(1:25,25,'Replace',true);
connICA_matrix=[connICA_matrix_control(sample_vector_control,:); connICA_matrix_epilepsy(sample_vector_epilepsy,:)];
%% 2. Setup configs 
addpath(genpath('FastICA_25'));  % add FastICA folder to path. Please check folder_name and update if needed.
configs.numRegions = 148; % aparc2009 parcellation
configs.mask = triu(true(configs.numRegions,configs.numRegions),1); % upper triangular are ones in 148x148 matrix
configs.numConn = size(connICA_matrix,1);%number of connectomes (subjects or sessions)
configs.numEdges = nnz(configs.mask);
configs.epsilon = 0.0001; configs.maxFinetune=100; configs.bootstrap=1;

configs.minFreq = 0.75;          % minimum frequency required 
configs.corrMin = 0.75;          % correlation between FC_traits
configs.corrMin_A = 0.75;        % correlation between subject weights, set it to 0 if bootstrapping subjects across numRuns 

%% ICA configs - update if necessary
configs.maxNumIterations=300; % number of times max. algorithm runs to find convergence
configs.numRuns=50;           % decreased from 100 to 50 to save computation time. Enrico used 100. I also used 100 for the original sample, no resampling.
configs.numRuns_max=200;      % number of runs maximum to obtain 50 converged runs

%% 3. Explore possible PCA dimensions -> grid dimensions
    % Each bootstrapped connICA has a different PCA structure, so recalculate
    %    how many eiganvectors we want to keep for the desired % of variance explained for each connICA
pca_var=[0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95];
numPCAComps_max = size(connICA_matrix,1);

% What could be the dimensions of the RC grid matrix? 
% Max number of eigenvalues for the 95% variance reduction - 15

numb_eigenvect_95PCA=zeros(1,100);
for j=1:100
    sample_vector_control=datasample(1:21,21,'Replace',true);
    sample_vector_epilepsy=datasample(1:25,25,'Replace',true);
    connICA_matrix=[connICA_matrix_control(sample_vector_control,:); connICA_matrix_epilepsy(sample_vector_epilepsy,:)];
    name=sprintf('bootstrap%d.mat',j);
    configs.which_dataset=name; 
    %configs.sample_vector=[sample_vector_control sample_vector_epilepsy];
    %configs.connICA_matrix=connICA_matrix;

    [~, ~, latent] = pca(connICA_matrix','NumComponents',numPCAComps_max); 
    variance = cumsum(latent)./sum(latent); 
    
    pca_eigenvectors=zeros(1,numel(pca_var));
    for i_pca=1:numel(pca_var)
        pca_eigenvectors(1,i_pca)=find(variance>pca_var(i_pca),1);
    end
    
    %RCs=cell(pca_eigenvectors(end),numel(pca_var)); % the max number of ICs (max number of PCA eigenvectors) = nrows in the RC grid
    
    numb_eigenvect_95PCA(j)=pca_eigenvectors(end);
  
end

% number of eigenvectors needed for 95% variance reduction in PCA, on the bootstrap sample with replacement
% varies between 19 and 30, based on how frequently the same person was
% selected in the with replacement resampling
disp(min(numb_eigenvect_95PCA)) % 19 
disp(max(numb_eigenvect_95PCA)) % 30

%% Resampling with replacements: Bootstrap
%% 4.a Parfor at the number of Independent components level (rows in the grid)
% keep in mind RCs change in size (hence the step before), typical size
% could be 26x12, only the number of rows change. Columns are the pca variances that are fix.
for j=1:1 %100 % specify number of times to obtain bootstrapped RCs table, ideally 100

    sample_vector_control=datasample(1:21,21,'Replace',true);
    sample_vector_epilepsy=datasample(1:25,25,'Replace',true);
    connICA_matrix=[connICA_matrix_control(sample_vector_control,:); connICA_matrix_epilepsy(sample_vector_epilepsy,:)];
    name=sprintf('bootstrap%d.mat',j);
    configs.which_dataset=name; 
    %configs.sample_vector=[sample_vector_control sample_vector_epilepsy];
    %configs.connICA_matrix=connICA_matrix;

    [~, ~, latent] = pca(connICA_matrix','NumComponents',numPCAComps_max); 
    variance = cumsum(latent)./sum(latent); 
    
    pca_eigenvectors=zeros(1,numel(pca_var));
    for i_pca=1:numel(pca_var)
        pca_eigenvectors(1,i_pca)=find(variance>pca_var(i_pca),1);
    end
 
    RCs=cell(pca_eigenvectors(end),numel(pca_eigenvectors)); % the max number of ICs (max number of PCA eigenvectors) = nrows in the RC grid
    
    tic
    for i=1:numel(pca_eigenvectors)
        lastEig=pca_eigenvectors(i);
        
        parfor k=2:lastEig
            [RC, RC_Index]=run_global_robust_connICA_parallel(connICA_matrix, configs, k, lastEig);
            RCs{k,i}=RC;  
        end
        %RCs{lastEig+1,i}=sprintf('elapsed time for this column %d',toc); %processor 35% Memory 65%
    end
    time_parfor_numOfIC=toc; %2.6501e+04
    filename=sprintf('bootstrap%d_RCs_replacement.mat',j);
    save(['C:\Users\user\Documents\Enrico\12_Consciousness_Connectomics\Bootstrap\' filename],'RCs');

end
clc; disp('done')
%% 4.b Parfor at both number of IC and PCA - create a linear grid vector(s) and use on for loop only
storage_grid_inf=struct;
for j=1:1%100 % specify number of times to resample and obtain the grid - ideally 100

    sample_vector_control=datasample(1:21,21,'Replace',true);
    sample_vector_epilepsy=datasample(1:25,25,'Replace',true);
    connICA_matrix=[connICA_matrix_control(sample_vector_control,:); connICA_matrix_epilepsy(sample_vector_epilepsy,:)];
    %name=sprintf('bootstrap%d.mat',j);
    %configs.which_dataset=name; 
    %configs.sample_vector=[sample_vector_control sample_vector_epilepsy];
    %configs.connICA_matrix=connICA_matrix;

    [~, ~, latent] = pca(connICA_matrix','NumComponents',numPCAComps_max); 
    variance = cumsum(latent)./sum(latent); 
    
    pca_eigenvectors=zeros(1,numel(pca_var));
    for i_pca=1:numel(pca_var)
        pca_eigenvectors(1,i_pca)=find(variance>pca_var(i_pca),1);
    end
    %
    grid1D=[]; % example
    grid1D(1,:)=repelem(pca_eigenvectors(1),(pca_eigenvectors(1)-1));  % -1 because we count from 2:pca_eigenvectors, don't calculate models with 1 independent component
    grid1D(2,:)=[2:pca_eigenvectors(1)];
    %
    for i=2:numel(pca_eigenvectors)
        grd=[]; % example
        grd(1,:)=repelem(pca_eigenvectors(i),(pca_eigenvectors(i)-1));  % -1 because we count from 2:pca_eigenvectors, don't calculate models with 1 independent component
        grd(2,:)=[2:pca_eigenvectors(i)];
        
        grid1D=[grid1D grd];
        
    end
    
    storage_grid_inf(j).grid1D=grid1D;
    storage_grid_inf(j).configs=configs;

    RCs=cell(1,size(grid1D,2)); % the max number of ICs (max number of PCA eigenvectors) = nrows in the RC grid
    %
    tic
    parfor i=1:size(grid1D,2)
        parameters=grid1D(:,i); 
        
        lastEig=parameters(1);
        numb_ic=parameters(2);
        [RC, RC_Index]=run_global_robust_connICA_parallel(connICA_matrix, configs, numb_ic, lastEig);
        RCs{i}=RC;
    end
    toc_grid1D=toc;
    disp(toc)    
    %filename=sprintf('bootstrap%d_RCs_other.mat',j);
    %save(['C:\Users\user\Documents\Enrico\12_Consciousness_Connectomics\Bootstrap\' filename],'RCs');

end
%clc; disp('done')
%% Notes
% 4.b: Running until 30 loops out of the 200 length entire loop. 100% processor
% use, 740 seconds computational time. Non-parallel simple for loop 777
% seconds.
% For the entire loop 36522 seconds, parfor.
% load bootstrap1_RCs.mat
% filename=sprintf('bootstrap%d_RCs_withrepl.mat',j);
% save(['C:\Users\user\Documents\Enrico\12_Consciousness_Connectomics\Bootstrap\' filename],'RCs');
%% 4.c regular for loop
for j=1:1 %100 % specify number of times to obtain bootstrapped RCs table, ideally 100

    sample_vector_control=datasample(1:21,21,'Replace',true);
    sample_vector_epilepsy=datasample(1:25,25,'Replace',true);
    connICA_matrix=[connICA_matrix_control(sample_vector_control,:); connICA_matrix_epilepsy(sample_vector_epilepsy,:)];
    name=sprintf('bootstrap%d.mat',j);
    configs.which_dataset=name; 
    %configs.sample_vector=[sample_vector_control sample_vector_epilepsy];
    %configs.connICA_matrix=connICA_matrix;

    [~, ~, latent] = pca(connICA_matrix','NumComponents',numPCAComps_max); 
    variance = cumsum(latent)./sum(latent); 
    
    pca_eigenvectors=zeros(1,numel(pca_var));
    for i_pca=1:numel(pca_var)
        pca_eigenvectors(1,i_pca)=find(variance>pca_var(i_pca),1);
    end
 
    RCs=cell(pca_eigenvectors(end),numel(pca_eigenvectors)); % the max number of ICs (max number of PCA eigenvectors) = nrows in the RC grid
    
    tic
    for i=1:1%numel(pca_eigenvectors)
        lastEig=pca_eigenvectors(i);
        for k=2:lastEig
            configs.lastEig=lastEig;
            configs.numOfIC=k;
            [RC, RC_Index]=run_global_robust_connICA(connICA_matrix, configs);
            RCs{k,i}=RC;
        end
        %RCs{lastEig+1,i}=sprintf('elapsed time for this column %d',toc); %processor 35% Memory 65%
    end
    time_forloop=toc; %28830 seconds for one grid
    filename=sprintf('bootstrap%d_RCs_replacement.mat',j);
    save(['C:\Users\user\Documents\Enrico\12_Consciousness_Connectomics\Bootstrap\' filename],'RCs');

end