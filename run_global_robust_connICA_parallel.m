%% Run the entire robust connICA algorithm with fastICA
% numberOfIC and lstEig is removed from the configs and need to be specified outside the configs so that the loop can run in paralell
function [RC, RC_Index, rc] = run_global_robust_connICA_parallel(connICA_matrix, configs, numberOfIC, lstEig)

if lstEig<numberOfIC 
    RC=[];
    RC_Index=[];
    rc=[];
    clc;
    warning('Number of PC is smaller than number of IC. fastICA did not run with this ICA setup')
else 

warning('on')
lastwarn(''); 
c=1;     % run count when algorithm converged
c_max=0; % run count for all iterations
clearvars A icasig
while c ~= configs.numRuns+1 % be careful might run for very long
    [icasig_onerun,A_onerun,~] = fastica(connICA_matrix,'approach','symm','numOfIC',numberOfIC,'verbose','off',...
        'epsilon',configs.epsilon,'maxNumIterations',configs.maxNumIterations,...
        'maxFinetune',configs.maxFinetune, 'lastEig', lstEig);%running fastica
    [warnMsg, warnId] = lastwarn;
    
    c_max=c_max+1;
    if c_max==configs.numRuns_max
        break
    end
    
    if isempty(warnMsg) % store those runs only where the algorithm converged
        A(:,:,c) = single(A_onerun);% connICA weights
        icasig(:,:,c) = single(icasig_onerun'); % connICA single run traits
        c=c+1;
    end
  
    lastwarn('');
    %if mod(c,25)==0
        %disp(sprintf('%d fastICA runs out of %d',c,configs.numRuns))
    %end
end

if c_max==configs.numRuns_max 
    RC=['Reached numRuns_max ', num2str(configs.numRuns_max), '. ', num2str(c) ,' successful runs out of desired ', num2str(configs.numRuns)];
    RC_Index=[];
    rc=['Reached numRuns_max ', num2str(configs.numRuns_max), '. ', num2str(c) ,' successful runs out of desired ', num2str(configs.numRuns)];
    clc;
    warning(['fastICA did not find results enough times with this ICA setup.', num2str(c) ,' successful runs out of desired ', num2str(configs.numRuns)])
else % Continue with everything if we succeeded with the previous step.

    
% Robust trait extraction:dr Enrico Amico 2017
[comp_match_run1, freq] = run_robust_connICA_fast_parallel(A,icasig,configs, numberOfIC, lstEig); % comp_mach_run1 stores the robust component ID per connICA run
aux=(freq>configs.minFreq); 

RC_Index = find(aux(:,1)); % Robust Components Index the numbers of components we keep.
 
if isempty(RC_Index)
    RC=['No robust traits'];
    RC_Index=[];
    rc=['No robust traits'];
    warning('No Robust FC traits found!')
    return
else
    RC = struct;
    rc = struct; % includes Ratio of Robust and Weights variance information only, for Comb. Ind.
for t=1:length(RC_Index) % 1 to 5, or as many as robust traits find previously
    compIndex=RC_Index(t); % choose the component(order of component) we want to look at 
    icasig_comp = nan(configs.numRuns,size(connICA_matrix,2));
    a0 = A(:,comp_match_run1(compIndex,1),1); %this is used as reference so that the same component across runs is always positively correlated (not inverted)
    weights = nan(size(A,1),configs.numRuns);
    for i=1:configs.numRuns 
        if comp_match_run1(compIndex,i)>0 % take only if robust index is larger than 1
            a = A(:,comp_match_run1(compIndex,i),i);
            icasig_comp_one = squeeze(icasig(:,comp_match_run1(compIndex,i),i));
            if corr(a,a0)>0
                icasig_comp(i,:) = icasig_comp_one;
            else
                a = -a;
                icasig_comp(i,:) = -icasig_comp_one;
            end
            weights(:,i) = a;  % weights for component t and run i are stored
        end
    end

    clear comp;
    % Calculations are commented because no need for this information, for the Combined Index at this point, uncomment if necessary
%     RC(t).vector = [];
%     RC(t).vector = nanmean(icasig_comp); %avg per column (across runs)
%     RC(t).matrix = zeros(configs.numRegions,configs.numRegions);
%     RC(t).matrix(configs.mask)=RC(t).vector; %fill upper triangular
%     RC(t).matrix = RC(t).matrix + (RC(t).matrix'); %symmetrize matrix
%     RC=rmfield(RC, 'vector'); % decrease storage space
%     RC(t).weights = weights; % store Robust component weights 
      RC(t).trait_weights_var=nanvar(weights');
      RC(t).trait_weights_var_real=nanvar(weights')./abs(nanmean(weights'));
    
    
%     RC(t).weights_mean= nanmean(weights'); % to be confirmed
%     RC(t).median_weights_mean=median(RC(t).weights_mean);
%     RC(t).weights_std=nanstd(weights');
%     RC(t).weights_var=nanvar(weights');
%     
%     RC(t).weights_std_real=nanstd(weights')./abs(nanmean(weights'));
%     RC(t).weights_var_real=nanvar(weights')./abs(nanmean(weights'));
%     
%     RC(t).weights_var_trait=mean([RC(t).weights_var]);
%     RC(t).inf_trait=round(sum(abs(RC(t).vector)),2);
    
end
%RC(1).configs=configs
%rc.weights_var=mean([RC.trait_weights_var]);
rc.ratio_of_robust=length(RC_Index)/numberOfIC;
rc.weights_var_real=mean([RC.trait_weights_var_real]); % over all traits and all subjects the mean relative variation of IC weights
end
end
end