% Advanced topics in network neuroscience EE-619
% Lecturer: Dr. Enrico Amico
% Seventh exercise class, code by EA, December 4 2020
% In this session we'll learn how to:
% - Compute connectome-based independent components (connICA, Amico et al., 2017)
% The code that follows can be freely downloaded from the CONNplexity lab website, https://engineering.purdue.edu/ConnplexityLab/publications

%set(0,'DefaultFigureWindowStyle','docked'); % to dock figures within command window
set(0,'DefaultFigureWindowStyle','normal'); % to go back to normal mode
clearvars; close all;


% batch example connICA (Amico et al., NeuroImage 2017)
% Enrico Amico & Joaquin Goni, Purdue University
% version 1.1 Apr 10th 2018
% Updates on this version:
% - Now saving the vector of weights in the RC structure
% - Displaying error (warning) messages when PCA comps retained is lower
% (equal) than ICA number of comps requested and set in configs.numOfIC
% PLEASE CITE US!
% If you are using this code for your research, please kindly cite us:
% Mapping the functional connectome traits of levels of consciousness.
% Amico, E. ... & Goñi, J. (2017). NeuroImage, 148, 201-211.
% http://www.sciencedirect.com/science/article/pii/S1053811917300204
%
% IMPORTANT: FastICA package IS NEEDED!
% Please download FastICA package
% https://research.ics.aalto.fi/ica/fastica/code/dlcode.shtml

%% initialize environment
close all;
clearvars
clc

%% Load FC data (functional connectivity). This is a sample data.
%  dimensions are brain regions x brain regions x subjects
load FC.mat; %164x164x20

%% Configuration
addpath(genpath('FastICA_25')); %add FastICA folder to path. Please check folder_name and update if needed.
% ConnICA params
configs.numRegions = 164; % aparc2009 parcellation. 164x164 FC.
configs.mask = triu(true(configs.numRegions,configs.numRegions),1); % upper triangular mask
configs.numConn = size(FC,3);% 20 pieces of FC matrices. Number of connectomes (subjects or sessions)
configs.numEdges = nnz(configs.mask); %13366 = 164*0.5*163

%% create connICA input matrix (nSubj x n functional edges)
connICA_matrix = zeros(configs.numConn,configs.numEdges); 
for i=1:configs.numConn % 20
    aux = FC(:,:,i);
    connICA_matrix(i,:) = aux(configs.mask); % Extract the FC as a vector of numbers
    % Each row is a vector of the upper triangular FC matrix corresponding to one subject
end
initial_connICA_matrix=connICA_matrix;
%% FastICA params

% PCA computes eigenvectors of the covariance matrix ("principal axes") 
%    and sorts them by their eigenvalues (amount of explained variance)

% ICA hyperparameters
configs.epsilon = 0.0001; % epsilon default is 0.0001
configs.maxNumIterations = 1000; % maxNumIterations default is 1000
configs.maxFinetune = 100; % maxFinetune default is 100
configs.numRuns = 100; %100; %run FastICA several time for robust comps extraction  
configs.numOfIC = 5; % number independent components (explore different range)    

flags.PCA = 1; % if 1, perform PCA before ICA
configs.PCAthr = 0.65; % Explained variance threshold for PCA. Feel free to explore 0.75 or .80
% perform PCA before ICA
if flags.PCA==1
    disp('running PCA compression before ICA...');
    numFCs = size(connICA_matrix,1);
    numPCAComps = size(connICA_matrix,1);
    [~, ~, latent] = pca(connICA_matrix','NumComponents',numPCAComps); 
    variance = cumsum(latent)./sum(latent); 
    variance = variance(1:numPCAComps); %explained variance with the selected num of PCA comps
    numPCAComps = find(variance>=configs.PCAthr,1);
    disp('# PCA comps retained');
    disp(numPCAComps);
    if numPCAComps<configs.numOfIC
        error('Error: PCA comps retained lower than number of requested ICA comps! Check the configs.numOfIC and try again');
    end
    if numPCAComps==configs.numOfIC
        warning(sprintf('Warning: number of PCA comps retained equal to the number of requested ICA comps.\n Please note that, in our experience, sometimes FastICA does not reach convergence under this configuration')); %#ok<SPWRN>
    end
    % Reduce to 9 principal components instead of 20 number of subjects
    % Obtain 13366x9 matrix - SCORE
    
    % 20x9 matrix - COEFFS
    %figure, plot(1:numFCs,variance,'ok'); % uncomment this line if you want to see the variance explained plot
    [COEFFS, SCORE, latent] = pca(connICA_matrix','NumComponents',numPCAComps); % When 'NumComponents' is specified pca gives the first k components only    
    % Latent gives back the eigenvalue of the covariance matrix of connICA
    % Keep in mind that PCA matlab command automatically centers the data
    %   by subtracting the column wise mean before the svd dec.
    PCA_clean_matrix = SCORE * COEFFS'; % PCA reconstructed de-meaned data
    
    PCA_clean_matrix = bsxfun(@plus, PCA_clean_matrix,mean(connICA_matrix,2)'); % Take the mean of each column
    % Elementwise addition to the matrix. 
    PCA_clean_matrix = PCA_clean_matrix'; %back to subjects x edges form. Now we have 9 instead of 20.
    connICA_matrix = PCA_clean_matrix; % overwrite the connICA with the PCA matrix
    disp('Done.');
end


%%
%save 'connICA_matrix_PCA9.mat' connICA_matrix
%%
icasig = nan(size(connICA_matrix,2),configs.numOfIC,configs.numRuns);
A = nan(size(connICA_matrix,1),configs.numOfIC,configs.numRuns);
for i=1:configs.numRuns
    [icasig_onerun,A_onerun,~] = fastica(connICA_matrix,'approach','symm','numOfIC',configs.numOfIC,'verbose','off',...
        'epsilon',configs.epsilon,'maxNumIterations',configs.maxNumIterations,...
        'maxFinetune',configs.maxFinetune);%running fastica
    A(:,:,i) = single(A_onerun);% connICA weights
    icasig(:,:,i) = single(icasig_onerun'); % connICA single run traits
    if mod(i,25)==0 % just to see how the loop goes, where we are
        disp(sprintf('%d runs out of %d',i,configs.numRuns)) %#ok<DSPS>
    end
end

%% robust traits extraction criteria (see Amico et al., 2017)
configs.minFreq = 0.75; % minimum frequency required 
configs.corrMin = 0.75; % correlation between FC_traits
configs.corrMin_A = 0.75; % correlation between subject weights, set it to 0 if bootstrapping subjects across numRuns 
[comp_match_run1, freq] = run_robust_connICA_fast(A,icasig,configs); % comp_mach_run1 stores the robust component ID per connICA run
aux=(freq>configs.minFreq);
RC_Index = find(aux(:,1)); % RC_Index == Robust Components Index
disp(RC_Index)
if isempty(RC_Index)
    error('No Robust FC traits found!');
end
%% Put back the robust traits in matrix form

load('aparc_a2009_yeo_RS7.mat');
RC = struct;
for t=1:length(RC_Index)
    compIndex=RC_Index(t); % choose the component(order of component) we want to look at 
    figure,
    icasig_comp = nan(configs.numRuns,size(connICA_matrix,2));
    a0 = A(:,comp_match_run1(compIndex,1),1); %this is used as reference so that the same component across runs is always positively correlated (not inverted)
    weights = nan(size(A,1),configs.numRuns);
    for i=1:configs.numRuns
        if comp_match_run1(compIndex,i)>0
            a = A(:,comp_match_run1(compIndex,i),i);
            icasig_comp_one = squeeze(icasig(:,comp_match_run1(compIndex,i),i));
            if corr(a,a0)>0
                plot(a); hold on;
                icasig_comp(i,:) = icasig_comp_one;
            else
                plot(-a); hold on;
                a = -a;
                icasig_comp(i,:) = -icasig_comp_one;
            end
            weights(:,i) = a;  % weights for component t and run i are stored
        end
    end
    ylabel('weights');
    xlabel('subjects');
    title(sprintf('connICA comp %d',compIndex));
    clear comp;
    RC(t).vector = [];
    RC(t).vector = nanmean(icasig_comp); %avg per column (across runs)
    RC(t).matrix = zeros(configs.numRegions,configs.numRegions);
    RC(t).matrix(configs.mask)=RC(t).vector; %fill upper triangular
    RC(t).matrix = RC(t).matrix + (RC(t).matrix'); %symmetrize matrix
    RC(t).weights = weights; % store Robust component weights per run
    
    figure; 
    imagesc(RC(t).matrix(yeoOrder,yeoOrder),[-3,3]); colormap jet; colorbar; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); xlabel('regions'); ylabel('regions');
    title(sprintf('connICA comp %d',compIndex));
end
