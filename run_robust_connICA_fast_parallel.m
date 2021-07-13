%% [comp_match_run1, freq] = run_robust_connICA_fast(A,icasig,configs)
% Extract the most robust ICA traits across runs, as explained in Amico et
% al., Neuroimage 2017.
% In summary, for the trait to be considered as robust it has to correlate higher than a certain threshold across runs both in terms of trait and weights (defined by user in configs.corrMin and configs.corrMin_A); 
% also it has to appear at least in a predefined of the total runs (set by user in configs.minFreq).
% Inputs:
% -A : 3D vector of ICA weights
% -icasig : 3D vector of ICA traits
% -configs: structure of params.
    % - numberOfIC: number of ICA components
    % - configs.numRuns: numnber of ICA runs
    % - configs.corrMin: number of minimum correlation between single-run vectors of traits for robustness selection  
    % - configs.corrMin_A = number of minimum correlation between single-run  vectors of weights for robustness selection
    % - configs.minFreq = number of minimum frequency across runs for robustness selection
% 
% -Outputs
% The robust traits are ordered and (if necessary) rearranged according to the first (best) run
% -comp_match_run1 = vector of the indices of the robust traits across runs 
% -freq = frequency of occurence of traits per runs
%
% example of use:
% [comp_match_run1, freq] = run_robust_connICA_fast(A,icasig,configs)
% Enrico Amico & Joaquin Goni, Purdue University
% version 1.0 Jan 16th 2018
function [comp_match_run1, freq] = run_robust_connICA_fast_parallel(A,icasig,configs, numberOfIC, lstEig)

% reshape into 2D vecs 
A_2D = single(reshape(A,size(A,1),numberOfIC*configs.numRuns)); % 21 , 50 The matrices are put next to each other from left to right
icasig_2D = single(reshape(icasig,size(icasig,1),numberOfIC*configs.numRuns));
% Compute correlation all vs all
%disp ('Computing corr btw comps (might take some time)..')
C_A = single(corr(A_2D,A_2D));
C_icasig = single(corr(icasig_2D,icasig_2D));
if (configs.corrMin<0.5) || ((configs.corrMin_A<0.5) && (configs.corrMin_A>0))
    warning(sprintf('Correlation threshold between traits set below 0.5.\n Please note that a higher threshold is recommended...'))
end
mask_freq = (abs(C_A) > configs.corrMin_A) & (abs(C_icasig) > configs.corrMin); % if weights AND icasig are > min Corr value

if (nnz(sum(mask_freq)>configs.numRuns)>0)
    error(sprintf('Error: trait co-occurence found more than once in at least one run!\n Please increase correlation thresholds:\n configs.corrMin and/or configs.corrMin_A'))
end
% reshape freq back to (numComp, numRuns)
freq = reshape(sum(mask_freq),numberOfIC,configs.numRuns);
freq = freq./configs.numRuns;


robust_comps_runs = sum(freq>configs.minFreq);
reference_run = find(robust_comps_runs==max(robust_comps_runs),1);
if reference_run ~=1 % change reference run
    Aaux = A(:,:,1);
    A(:,:,1) = A(:,:,reference_run);
    A(:,:,reference_run) = Aaux;
    icasig_aux = icasig(:,:,1);
    icasig(:,:,1) = icasig(:,:,reference_run);
    icasig(:,:,reference_run) = icasig_aux;
    
    % Update corr mats
    %disp ('New ref run: updating (might take some time)..')
    A_2D = single(reshape(A,size(A,1),numberOfIC*configs.numRuns));
    icasig_2D = single(reshape(icasig,size(icasig,1),numberOfIC*configs.numRuns));
    C_A = single(corr(A_2D,A_2D));
    C_icasig = single(corr(icasig_2D,icasig_2D));
    % update freq mat
    mask_freq = (abs(C_A) > configs.corrMin_A) & (abs(C_icasig) > configs.corrMin); % if weights AND icasig are > min Corr value
    % reshape freq back to (numComp, numRuns)
    freq = reshape(sum(mask_freq),numberOfIC,configs.numRuns);
    freq = freq./configs.numRuns;
    % compute comp_match_run1
    comp_match_run1 = zeros(numberOfIC,configs.numRuns);
    comp_match_run1(:,1) = 1:numberOfIC;
    for i=1:numberOfIC
        for k=2:configs.numRuns
            range = ((k-1)*numberOfIC)+1:((k-1)*numberOfIC)+numberOfIC;
            c = abs(C_A(i,range));
            ctrait = abs(C_icasig(i,range));
            if max(c) > configs.corrMin_A && max(ctrait) > configs.corrMin
                pos = find(ctrait==max(ctrait));
                comp_match_run1(i,k) = pos;
            end
        end
    end    
else
    % compute comp_match_run1
    comp_match_run1 = zeros(numberOfIC,configs.numRuns);
    comp_match_run1(:,1) = 1:numberOfIC;
    for i=1:numberOfIC
        for k=2:configs.numRuns
            range = ((k-1)*numberOfIC)+1:((k-1)*numberOfIC)+numberOfIC;
            c = abs(C_A(i,range));
            ctrait = abs(C_icasig(i,range));
            if max(c) > configs.corrMin_A && max(ctrait) > configs.corrMin
                pos = find(ctrait==max(ctrait));
                comp_match_run1(i,k) = pos;
            end
        end
    end    
end
