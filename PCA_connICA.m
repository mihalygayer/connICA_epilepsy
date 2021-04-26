%% Set up the connICA in another file
connICA_matrix=connICA_matrix_epilepsy;
%% Running PCA 
flags.PCA=1
configs.PCAthr = 0.9; 
% perform PCA before ICA
if flags.PCA==1
    numFCs = size(connICA_matrix,1);
    numPCAComps = size(connICA_matrix,1);
    [~, ~, latent] = pca(connICA_matrix','NumComponents',numPCAComps); 
    variance = cumsum(latent)./sum(latent); 
    variance = variance(1:numPCAComps); %explained variance with the selected num of PCA comps
    numPCAComps = find(variance>=configs.PCAthr,1);
    disp('# PCA comps retained');
    disp(numPCAComps);
    close all;
    
    %%%%%%%%%%%%%%%%%%%%%  variance explained plot, update title name
    figure, plot(1:numFCs,variance,'ok'); xlabel('number of principal components'); ylabel('Variance explained');
    title('Variance explained by number of PC: sample dataset')
    
    [COEFFS, SCORE, latent] = pca(connICA_matrix','NumComponents',numPCAComps); % When 'NumComponents' is specified pca gives the first k components only    
    % Latent gives back the eigenvalue of the covariance matrix of connICA
    % Keep in mind that PCA matlab command automatically centers the data
    %   by subtracting the column wise mean before the svd dec.
    PCA_clean_matrix = SCORE * COEFFS'; % PCA reconstructed de-meaned data
    
    PCA_clean_matrix = bsxfun(@plus, PCA_clean_matrix,mean(connICA_matrix,2)'); % Take the mean of each column
    % Elementwise addition to the matrix. 
    PCA_clean_matrix = PCA_clean_matrix'; %back to subjects x edges form. Now we have 9 instead of 20.
    %connICA_matrix = PCA_clean_matrix; % overwrite the connICA with the PCA matrix
    disp('Done.');
end
%%

%%
pca_result=vertcat(round(variance,3)', round(1:numFCs,1))
%%
save 'pca_result_sample.mat' pca_result