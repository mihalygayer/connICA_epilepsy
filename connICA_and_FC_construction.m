% Master Project, Mihály Gayer
% Construct the connICA matrix and FCs for healthy and the epileptic patients. 
% Data from Dr. Jonathan Wirsich, Slack, 04.2021, and Dr. Enrico Amico (Healthy control class dataset)
% 1. Import: - healthy control dataset (Zenodo, Dr. Wirshich, 256Ch3T)
%            - the FC from the healthy control (Class dataset from Dr. Amico) 
%            - the two epileptic datasets, rtle and ltle separately
%            
%    Apply for-loop to extract data subject by subject to construct the FC, 
%    from the FC construct the connICA matrix

% 2. Merge and save different combination of connICA matrices

% 3. Plot the FC result using imagesc 
% 4. Obtain groupwise mean functional connectome for each 4 groups
%%
clearvars;
clc; addpath('C:\Users\user\Documents\Enrico\12_Consciousness_Connectomics')

%% Configuration - Dr. Enrico Amico code from Class of N.N. 2020
configs.numRegions=148; % There are 148 regions in this parcellation
configs.mask = triu(true(configs.numRegions,configs.numRegions),1); % to extract upper triangular part of the matrix for connICA

%% Healthy control Zenodo, 256Ch3T - 21 subjects
% Healthy control data is downloaded from: https://zenodo.org/record/3905103#.YBgYtOhKg2z
load('eeg-fmri_connectomes_destrieux_scrubbed_256Ch3T.mat')

% aim: have the same data structure as the one in Dr. Amico's code for FC matrix -> connICA matrix
FC_control=struct; % Construct the FC matrices for each subject - Dr. Jonathan Wirshich
for i = 1:length(subj)
 mrtx = zeros(148); 
 count = 1;
 for r1 = 1:148-1
  for r2 = r1+1:148
   mrtx(r1, r2) = subj(i).sess.fMRI(count);
   mrtx(r2, r1) = mrtx(r1, r2);
   count = count + 1;
  end
 end
 FC_control(i).matrix=mrtx;
end

connICA_matrix_control=zeros(length(subj),length(subj(1).sess.fMRI));
for i = 1:length(subj)
    aux = FC_control(i).matrix(:,:);
    connICA_matrix_control(i,:) = aux(configs.mask); % Extract the FC matrix as a vector, using upper triangular of ones
end
%% save/load
save 'connICA_matrix_control.mat' connICA_matrix_control
save 'FC_control.mat' FC_control

%% Healthy control CLASS dataset from Dr. Amico N.N 2020
% FC Class dataset - Advanced topics in network neuroscience 2020, EPFL E-619
load FC.mat % 164x164x148

% Transfrom FC from 164x164 -> 148x148
FC_class=struct;
for i=1:size(FC,3)
    FC_class(i).matrix=FC([1:74 90:163],[1:74 90:163],i); % 2021.02.01 Enrico's comment on which to remove
    FC_class(i).matrix=FC_class(i).matrix-diag(diag(FC_class(i).matrix)); % make diagonal 0
end

connICA_matrix_class = zeros(size(FC,3),length(subj(1).sess.fMRI)); 
for i=1:size(FC,3)
    aux = FC_class(i).matrix;
    connICA_matrix_class(i,:) = aux(configs.mask); 
end
%%
save 'connICA_matrix_class.mat' connICA_matrix_class
save 'FC_class.mat' FC_class

%% rtle - 11 subjects
% Dr. Jonathan Wirsich, 2021.03
load eeg-fmri_connectomes_destrieux_scrubbed_truncTo5min_256Ch3T_eeg-fmri-geneva_rtle_notimeslice.mat

FC_rtle=struct;
for i = 1:length(subj)
 mrtx = zeros(148); 
 count = 1;
 for r1 = 1:148-1
  for r2 = r1+1:148
   mrtx(r1, r2) = subj(i).sess.fMRI(count);
   mrtx(r2, r1) = mrtx(r1, r2);
   count = count + 1;
  end
 end
 FC_rtle(i).matrix=mrtx;
end

% Amico
connICA_matrix_epilepsy_rtle=zeros(length(subj),length(subj(1).sess.fMRI));
for i=1:length(subj)
    aux = FC_rtle(i).matrix(:,:);
    connICA_matrix_epilepsy_rtle(i,:)=aux(configs.mask);
end

%% llte - 14 subjects
% Dr. Jonathan Wirsich, 2021.03
load eeg-fmri_connectomes_destrieux_scrubbed_truncTo5min_256Ch3T_eeg-fmri-geneva_ltle_notimeslice.mat

FC_ltle=struct;
for i = 1:length(subj)
 mrtx = zeros(148); 
 count = 1;
 for r1 = 1:148-1
  for r2 = r1+1:148
   mrtx(r1, r2) = subj(i).sess.fMRI(count);
   mrtx(r2, r1) = mrtx(r1, r2);
   count = count + 1;
  end
 end
 FC_ltle(i).matrix=mrtx;
end

connICA_matrix_epilepsy_ltle=zeros(length(subj),length(subj(1).sess.fMRI));
for i=1:length(subj)
    aux = FC_ltle(i).matrix(:,:);
    connICA_matrix_epilepsy_ltle(i,:)=aux(configs.mask);
end

%% Both epilepsy, merge
connICA_matrix_epilepsy=[connICA_matrix_epilepsy_rtle; connICA_matrix_epilepsy_ltle];
FC=[FC_rtle, FC_ltle];
save 'connICA_matrix_epilepsy.mat' connICA_matrix_epilepsy
save 'FC_epilepsy.mat' FC

%% 2 healthy controls as one
connICA_matrix_2healthy=[connICA_matrix_control; connICA_matrix_class];
save 'connICA_matrix_2healthy.mat' connICA_matrix_2healthy

%% Sample = Jonathan control + Epilepsy 
load connICA_matrix_control.mat % Jonathan control dataset
connICA_matrix_sample=[connICA_matrix_control; connICA_matrix_epilepsy];
save 'connICA_matrix_sample.mat' connICA_matrix_sample

%% Class healthy + epilepsy
connICA_matrix_class_epilepsy=[connICA_matrix_class; connICA_matrix_epilepsy];
save 'connICA_matrix_class_epilepsy.mat' connICA_matrix_class_epilepsy
%% All = Class + Jonathan control + Epilepsy
load connICA_matrix_class.mat
load connICA_matrix_control.mat 
%%
connICA_matrix_all=[connICA_matrix_class; connICA_matrix_control; connICA_matrix_epilepsy];
save 'connICA_matrix_all.mat' connICA_matrix_all

%%%%% 3.
%% Display Functional Connectomes of epileptic patients and healthy controls
load('aparc_a2009_yeoRS7_148reg_eeg_nosubc.mat'); % from Jonathan's email, 2021.02.01
%% rtle - Display FC in loop, 3 rows 4 columns
hold off
t=tiledlayout(3,4)
for i=1:length(FC_rtle)-1
    nexttile
    imagesc(FC_rtle(i).matrix(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title(sprintf('Patient %i',i), 'fontsize',8);
end
nexttile
    imagesc(FC_rtle(11).matrix(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title(sprintf('Patient %i',i+1), 'fontsize',8); colorbar;
title(t,'Functional Connectomes of right temporal lobe epileptic patients')
%% ltle - Display FC in loop, 3 rows 4 columns
hold off
t=tiledlayout(3,5)
for i=1:length(FC_ltle)-1
    nexttile
    imagesc(FC_ltle(i).matrix(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title(sprintf('Patient %i',i), 'fontsize',8);
end
nexttile
    imagesc(FC_ltle(14).matrix(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title(sprintf('Patient %i',i+1), 'fontsize',8); colorbar;
title(t,'Functional Connectomes of left temporal lobe epileptic patients')

%%
load FC_control
%%
hold off
t=tiledlayout(3,7)
for i=1:length(FC_control)-1
    nexttile
    imagesc(FC_control(i).matrix(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title(sprintf('Subject %i',i), 'fontsize',8);
end
nexttile
    imagesc(FC_control(21).matrix(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title(sprintf('Patient %i',i+1), 'fontsize',8); colorbar;
title(t,'Functional Connectomes of healthy control subjects (Zenodo)')

%% FC Class dataset - Advanced topics in network neuroscience 2020, EPFL E-619
load FC_class.mat
hold off
t=tiledlayout(3,7)
for i=1:length(FC_class)-1
    nexttile
    imagesc(FC_class(i).matrix(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title(sprintf('Subject %i',i), 'fontsize',8);
end
nexttile
    imagesc(FC_class(20).matrix(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title(sprintf('Subject %i',i+1), 'fontsize',8); colorbar;
title(t,'Functional Connectomes of healthy control subjects (Class)');

%%%% Mean functional connectomes by 4 groups
%% Mean FC_class
clc;
FC_class_mean=zeros(148,148);
for i=1:length(FC_class)
FC_class_mean=FC_class_mean+FC_class(i).matrix;
end
FC_class_mean=FC_class_mean/length(FC_class);

imagesc(FC_class_mean(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title('Average Functional connectome, healthy controls (class)', 'fontsize',10); colorbar;
%% Mean FC_control
FC_control_mean=zeros(148,148);
for i=1:length(FC_control)
    FC_control_mean=FC_control_mean+FC_control(i).matrix;
end
FC_control_mean=FC_control_mean/length(FC_control);
imagesc(FC_control_mean(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title('Average Functional connectome, healthy controls (Zenodo)', 'fontsize',10); colorbar;
    
%% Mean FC_rtle
FC_rtle_mean=zeros(148,148);
for i=1:length(FC_rtle)
    FC_rtle_mean=FC_rtle_mean+FC_rtle(i).matrix;
end
FC_rtle_mean=FC_rtle_mean/length(FC_rtle);
imagesc(FC_rtle_mean(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title('Average Functional connectome, epileptic patients (RTLE)', 'fontsize',10); colorbar;
    
%% Mean FC_rtle
FC_ltle_mean=zeros(148,148);
for i=1:length(FC_ltle)
    FC_ltle_mean=FC_ltle_mean+FC_ltle(i).matrix;
end
FC_ltle_mean=FC_ltle_mean/length(FC_ltle);
imagesc(FC_ltle_mean(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title('Average Functional connectome, epileptic patients (LTLE)', 'fontsize',10); colorbar;
    
%% Average functional connectome by group
a=tiledlayout(1,4)
nexttile
imagesc(FC_class_mean(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title('Healthy controls (class)', 'fontsize',8); %colorbar;
nexttile
imagesc(FC_control_mean(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title('Healthy controls (Zenodo)', 'fontsize',8); %colorbar;
nexttile
imagesc(FC_rtle_mean(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title('Epileptic patients (RTLE)', 'fontsize',8); %colorbar;
nexttile
imagesc(FC_ltle_mean(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title('Epileptic patients (LTLE)', 'fontsize',8); %colorbar;
title(a, 'Average functional connectome of each subject group')