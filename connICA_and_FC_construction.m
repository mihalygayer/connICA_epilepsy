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
% 5. Attach subject information imported from 3 tsv files one for control two for epileptic patients, as of 05.19.2021 from Dr.Wirshich
% 6. Explore participant information in each group
% 7. Extract FDA movement information for each subject
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
 FC_control(i).fd_power_mean=mean(subj(i).sess.fd_power); % Save 150 fd_movement parameter for each subject
end

connICA_matrix_control=zeros(length(subj),length(subj(1).sess.fMRI));
for i = 1:length(subj)
    aux = FC_control(i).matrix(:,:);
    connICA_matrix_control(i,:) = aux(configs.mask); % Extract the FC matrix as a vector, using upper triangular of ones
end

%% save/load
% save 'connICA_matrix_control.mat' connICA_matrix_control
% save 'FC_control.mat' FC_control

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
% save 'connICA_matrix_class.mat' connICA_matrix_class
% save 'FC_class.mat' FC_class

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
 FC_rtle(i).fd_movement_mean=mean(subj(i).sess.fd_power); % mean of 300 entries fd_movement parameter
end

% Amico's data transformation
connICA_matrix_epilepsy_rtle=zeros(length(subj),length(subj(1).sess.fMRI));
for i=1:length(subj)
    aux = FC_rtle(i).matrix(:,:);
    connICA_matrix_epilepsy_rtle(i,:)=aux(configs.mask);
end

%% llte - 14 subjects
% Dr. Jonathan Wirsich, 2021.03
load eeg-fmri_connectomes_destrieux_scrubbed_truncTo5min_256Ch3T_eeg-fmri-geneva_ltle_notimeslice.mat
%%
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
 FC_ltle(i).fd_movement_mean=mean(subj(i).sess.fd_power); % mean of 300 entries fd_movement parameter
end

% Amico's data transformation
connICA_matrix_epilepsy_ltle=zeros(length(subj),length(subj(1).sess.fMRI));
for i=1:length(subj)
    aux = FC_ltle(i).matrix(:,:);
    connICA_matrix_epilepsy_ltle(i,:)=aux(configs.mask);
end
%% Save 
% save FC_rtle.mat FC_rtle
% save FC_ltle.mat FC_ltle
%% Both epilepsy, merge
connICA_matrix_epilepsy=[connICA_matrix_epilepsy_rtle; connICA_matrix_epilepsy_ltle];
FC_epilepsy=[FC_rtle, FC_ltle];
% save 'connICA_matrix_epilepsy.mat' connICA_matrix_epilepsy
% save 'FC_epilepsy.mat' FC_epilepsy

%% 2 healthy controls as one
connICA_matrix_2healthy=[connICA_matrix_control; connICA_matrix_class];
% save 'connICA_matrix_2healthy.mat' connICA_matrix_2healthy

%% Sample = Jonathan control + Epilepsy 
load connICA_matrix_control.mat % Jonathan control dataset
connICA_matrix_sample=[connICA_matrix_control; connICA_matrix_epilepsy];
% save 'connICA_matrix_sample.mat' connICA_matrix_sample

%% Class healthy + epilepsy
connICA_matrix_class_epilepsy=[connICA_matrix_class; connICA_matrix_epilepsy];
% save 'connICA_matrix_class_epilepsy.mat' connICA_matrix_class_epilepsy
%% All = Class + Jonathan control + Epilepsy
load connICA_matrix_class.mat
load connICA_matrix_control.mat 
connICA_matrix_all=[connICA_matrix_class; connICA_matrix_control; connICA_matrix_epilepsy];
% save 'connICA_matrix_all.mat' connICA_matrix_all

%%%%% 3.
%% Display Functional Connectomes of epileptic patients and healthy controls
load FC_control.mat
load FC_rtle.mat
load FC_ltle.mat
load FC_epilepsy.mat

load('aparc_a2009_yeoRS7_148reg_eeg_nosubc.mat'); % from Jonathan's email, 2021.02.01
%% rtle - Display FC in loop, 3 rows 4 columns
hold off; clc;
close all
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
    title(sprintf('Patient %i',i+1), 'fontsize',8); 
    clrbar=colorbar; set(clrbar.XLabel,{'String','Rotation','Position'},{'Pearson correlation',0,[0.5 -1.2]})
    set(gcf,'color','w');
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
    title(sprintf('Patient %i',i+1), 'fontsize',8);
    clrbar=colorbar; set(clrbar.XLabel,{'String','Rotation','Position'},{'Pearson correlation',0,[0.5 -1.2]})
    set(gcf,'color','w');
title(t,'Functional Connectomes of left temporal lobe epileptic patients')

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
    title(sprintf('Patient %i',i+1), 'fontsize',8); 
    clrbar=colorbar; set(clrbar.XLabel,{'String','Rotation','Position'},{'Pearson correlation',0,[0.5 -1.2]})
    set(gcf,'color','w');
title(t,'Functional Connectomes of healthy control subjects')

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
    title(sprintf('Subject %i',i+1), 'fontsize',8);
    clrbar=colorbar; set(clrbar.XLabel,{'String','Rotation','Position'},{'Pearson correlation',0,[0.5 -1.2]})
    set(gcf,'color','w');
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
    
% Mean FC_rtle
FC_rtle_mean=zeros(148,148);
for i=1:length(FC_rtle)
    FC_rtle_mean=FC_rtle_mean+FC_rtle(i).matrix;
end
FC_rtle_mean=FC_rtle_mean/length(FC_rtle);
imagesc(FC_rtle_mean(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title('Average Functional connectome, epileptic patients (RTLE)', 'fontsize',10); colorbar;
    
% Mean FC_rtle
FC_ltle_mean=zeros(148,148);
for i=1:length(FC_ltle)
    FC_ltle_mean=FC_ltle_mean+FC_ltle(i).matrix;
end
FC_ltle_mean=FC_ltle_mean/length(FC_ltle);
imagesc(FC_ltle_mean(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title('Average Functional connectome, epileptic patients (LTLE)', 'fontsize',10); colorbar;
    
%% Average functional connectome by group
close all; a=tiledlayout(1,3)
% nexttile
% imagesc(FC_class_mean(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
%     set(gca,'xtick',[]); set(gca,'ytick',[]); 
%     title('Healthy controls (class)', 'fontsize',8); %colorbar;
nexttile
imagesc(FC_control_mean(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title('Healthy controls', 'fontsize',8); %colorbar;
nexttile
imagesc(FC_rtle_mean(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title('Epileptic patients (RTLE)', 'fontsize',8); %colorbar;
nexttile
imagesc(FC_ltle_mean(yeoOrder,yeoOrder),[-1,1]); colormap jet; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title('Epileptic patients (LTLE)', 'fontsize',8); %colorbar;
%title(a, 'Average functional connectome of each subject group');
    clrbar=colorbar; 
    set(clrbar.XLabel,{'String','Rotation','Position'},{'Pearson correlation',0,[0.5 -1.1]})
    set(gcf,'color','w');


%% 5. Subject information
% import manually the tsv files provided by Dr. Wirsich, then save it as a .mat file
%rTLEparticipants=rTLEparticpants;

%save 'rTLEparticipants.mat' rTLEparticipants
%save 'controlparticipants.mat' controlparticipants

% lTLEparticipants=lTLEparticpants;
% save lTLE_particpants.mat lTLEparticipants
%%
controlparticipants.participant_id_new=[1:21]';
controlparticipants=controlparticipants(:,[5 1 2 3 4]);
%%
rTLEparticipants.participant_id_new=[1:11]';
%%
rTLEparticipants=rTLEparticipants(:,[7 1 2 3 4 5 6]);
%%
lTLEparticipants.participant_id_new=[1:14]';
lTLEparticipants=lTLEparticipants(:,[7 1 2 3 4 5 6]);

%% 
epilepsy_participants=[rTLEparticipants; lTLEparticipants];
%save epilepsy_participants.mat epilepsy_participants
%% Participants information
clc; clearvars;
load controlparticipants.mat
load rTLE_participants.mat
load lTLE_particpants.mat
load epilepsy_participants.mat
%%

%% gender
clc;
sum(rTLEparticipants.sex=='M')/11
sum(lTLEparticipants.sex=='M')/14
(sum(rTLEparticipants.sex=='M')+sum(lTLEparticipants.sex=='M'))/25
sum(controlparticipants.sex=='M')/21
%% RTLE - age
clc;
min(rTLEparticipants.age)
max(rTLEparticipants.age)
mean(rTLEparticipants.age)
median(rTLEparticipants.age)
std(rTLEparticipants.age)
scatter(1:11, rTLEparticipants.age)
%% LTLE - age
clc;
min(lTLEparticipants.age)
max(lTLEparticipants.age)
mean(lTLEparticipants.age)
median(lTLEparticipants.age)
std(lTLEparticipants.age)

%% control - age
clc;
min(controlparticipants.age)
max(controlparticipants.age)
mean(controlparticipants.age)
median(controlparticipants.age)
std(controlparticipants.age)
%scatter(1:21, controlparticipants.age)

%% Epilepsy onset chart

close all;
t=tiledlayout(1,2)
nexttile
bar(sort(rTLEparticipants.epilepsy_onset), 'r'); xlabel('subject no.'); ylabel('epilepsy duration (years)')
title('RTLE patients')
nexttile
bar(sort(lTLEparticipants.epilepsy_onset), 'r'); xlabel('subject no.'); ylabel('epilepsy duration (years)')
title('LTLE patients')
title(t,'Epilepsy onset for each patient in the two epileptic groups')
set(gcf,'color','w');
%% RTLE - epilepsy onset
clc; 
min(rTLEparticipants.epilepsy_onset)
max(rTLEparticipants.epilepsy_onset)
mean(rTLEparticipants.epilepsy_onset)
median(rTLEparticipants.epilepsy_onset)
std(rTLEparticipants.epilepsy_onset)

%% LTLE - epilepsy onset
clc; 
min(lTLEparticipants.epilepsy_onset)
max(lTLEparticipants.epilepsy_onset)
mean(lTLEparticipants.epilepsy_onset)
median(lTLEparticipants.epilepsy_onset)
std(lTLEparticipants.epilepsy_onset)

%% Epilepsy onset both groups
mean([lTLEparticipants.epilepsy_onset; rTLEparticipants.epilepsy_onset])
std([lTLEparticipants.epilepsy_onset; rTLEparticipants.epilepsy_onset])
%% Epilepsy duration chart

close all;
t=tiledlayout(1,2)
nexttile
bar(sort(rTLEparticipants.epilepsy_duration), 'r'); xlabel('subject no.'); ylabel('epilepsy duration (years)')
title('RTLE patients')
nexttile
bar(sort(lTLEparticipants.epilepsy_duration), 'r'); xlabel('subject no.'); ylabel('epilepsy duration (years)')
title('LTLE patients')
title(t,'Epilepsy duration for each patient in the two epileptic groups')
set(gcf,'color','w');
%% RTLE - epilepsy duration
clc; 
min(rTLEparticipants.epilepsy_duration)
max(rTLEparticipants.epilepsy_duration)
mean(rTLEparticipants.epilepsy_duration)
median(rTLEparticipants.epilepsy_duration)
std(rTLEparticipants.epilepsy_duration)

%% LTLE - epilepsy duration
clc; 
min(lTLEparticipants.epilepsy_duration)
max(lTLEparticipants.epilepsy_duration)
mean(lTLEparticipants.epilepsy_duration)
median(lTLEparticipants.epilepsy_duration)
std(lTLEparticipants.epilepsy_duration)

