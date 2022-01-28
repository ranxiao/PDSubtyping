%% A standardized procedure for k-means clustering for PD subtyping and characterization of related attributes
% Complete data analysis codes for the study "Characterization of Parkinson’s Disease Subtypes and Related Attributes" by Shakya et al. 2022
% Author: Ran Xiao, Duke Univeristy, 2022
%% Part I: Baseline motor and non-motor data preparation
% load baseline data in PPMI
Data_bl = readtable('PPMI_all/PPMI_Baseline_Data_02Jul2018.csv');
% Select the following 20 variables
%   [Data_bl.PATNO Data_bl.APPRDX Data_bl.gen Data_bl.ageonset Data_bl.updrs3_score Data_bl.upsit Data_bl.moca Data_bl.gds Data_bl.stai Data_bl.NP1APAT Data_bl.NP1HALL Data_bl.NP1FATG ...
%    Data_bl.rem Data_bl.scopa Data_bl.scopa_gi Data_bl.scopa_ur Data_bl.scopa_cv Data_bl.scopa_therm Data_bl.scopa_pm Data_bl.scopa_sex];
Data_bl_onFile = Data_bl(:,[2 3 7 15 31 44 84 48 82 36 39 41 66 69:75]);
% select only PD patient data
Data_bl_onFile_pd = Data_bl_onFile(Data_bl_onFile.APPRDX==1,[1 3:end]);

% pain score, a non-motor symptom data used in the study, can be retrieved in the file below. 
UPDRS_part1 = readtable('PPMI_all/MDS_UPDRS_Part_I__Patient_Questionnaire'); % retreive pain score
% select 3 columns, including patNO, EventID (events include baseline,subsequent visits etc.), and pain score
pain_sc = UPDRS_part1(:,[3 4 10]);
% keep only rows with BL as the event_ID
pain_sc_bl = pain_sc(strcmp(pain_sc.EVENT_ID,'BL'),:);

% retrieving motor subscores in the following file
UPDRS_part3  = readtable('PPMI_all/MDS_UPDRS_Part_III'); 
% select items in UPDRS III to calculate UPDRS 4 motor subscores, 
% including tremor, regidity, bradykinesia and axial
motor_sc = UPDRS_part3(:,[3 4 12:47]); 
% derive tremor score
motor_sc.tremor = mean(motor_sc{:,[26:35]},2); 
% derive regidity score
motor_sc.regidity = mean(motor_sc{:,5:9},2);
%Bradykinesia Part III 2, 4-9, 14;
motor_sc.bradykinesia=mean(motor_sc{:, [4 10:20 25]}, 2);
%Axial subscore Part III 1, 9-13;
%In Li article both axial and bradykinesia have item 9 in common;
motor_sc.axial=mean(motor_sc{:, [3 20:24]}, 2);
% keep only rows with BL as event_ID
motor_sc_bl = motor_sc(strcmp(motor_sc.EVENT_ID,'BL'),[1 39:42]);

% adding pain score and 4 motor subscores to the "big" baseline data table
Baseline_allVar = innerjoin(Data_bl_onFile_pd,pain_sc_bl,'Keys','PATNO','RightVariables',3);
Baseline_allVar = innerjoin(Baseline_allVar,motor_sc_bl,'Keys','PATNO','RightVariables',2:5);

% remove patient data with missing values
tb_sum = sum(Baseline_allVar{:,2:end},2);
% find nan and remove corresponding patient data
Baseline_allVar_rm = Baseline_allVar(~isnan(tb_sum),:); 
% remove the first column representing patient ID, as the final data matrix
% for cluster analysis
dataBaseline = Baseline_allVar_rm(:,2:end);

% save data matrix
if ~exist('./results')
    mkdir('./results');
end
save('./results/dataBaseline', 'dataBaseline','Baseline_allVar_rm');
%% Part II: preparing the progression data
% load baseline data that contain symptom data at baseline visit, i.e. y0
Data_bl = readtable('PPMI_all/PPMI_Baseline_Data_02Jul2018.csv');
% for baseline, select only patient number, h&y scores at off state and diagnosis info
Data_bl_hy = [Data_bl.PATNO Data_bl.NHY Data_bl.updrs3_score Data_bl.moca Data_bl.APPRDX];
% select only PD patient data
Data_bl_hy = Data_bl_hy(Data_bl_hy(:,end)==1,1:end-1);

% load data from follow-up visits at year 1 to year 3
Data_y13 = readtable('PPMI_all/PPMI_Year_1-3_Data_02Jul2018.csv');
% for y1-3 follow-up, select  patient number, year of follow-up,  h&y scores at off state and diagnosis info
Data_y13_hy = [Data_y13.PATNO Data_y13.YEAR Data_y13.NHY Data_y13.updrs3_score Data_y13.moca Data_y13.APPRDX];
% select only PD patient data
Data_y13_hy = Data_y13_hy(Data_y13_hy(:,end)==1,1:end-1);

% this block of codes reorganize data into [patient number, y1_h&y,y2_h&y, y3_h&y]  
UniPat_y13 = unique(Data_y13_hy(:,1)); % find unique patients in the group
% initializing a variable with desired size, i.e., number of patients by 4 variables
Data_y13_hy_reorg = nan(length(UniPat_y13),10);
% loop through each patient data
for j =1: length(UniPat_y13) 
    %find all rows belong to each patient
    patData = Data_y13_hy(Data_y13_hy(:,1)==UniPat_y13(j),:); 
    % assign patient number to the first column
    Data_y13_hy_reorg(j,1) = UniPat_y13(j);
    % go through each row of patient-specific data and put y1_h&y to 2nd
    % column, y2_h&y to 3rd column, y3_h&y to 4th column
    for k = 1:size(patData,1)
        ind = patData(k,2);
        Data_y13_hy_reorg(j,ind+1) = patData(k,3);
        Data_y13_hy_reorg(j,ind+4) = patData(k,4);
        Data_y13_hy_reorg(j,ind+7) = patData(k,5);        
    end
end

% find unique patients in baseline data  
UniPat_bl = unique(Data_bl_hy(:,1));
% remove patients in baseline but without y1_3 h&y scores
Data_bl_hy_rm = Data_bl_hy(ismember(UniPat_bl,UniPat_y13),:);

% join baseline and y1_y3 data together
Data_joint_hy = [Data_bl_hy_rm Data_y13_hy_reorg];
% remove redundant patient-number column during the joining process
Data_joint_hy(:,5)=[];
% reorder the variable into the following: [patno, H&Y year0-3, UPDRS-III year0-3, and MOCA year0-3]
dataProgression = Data_joint_hy(:,[1 2 5 6 7 3 8 9 10 4 11 12 13]);

% save to results folder
save('./results/dataProgression', 'dataProgression');
%% Part III: a standardized and automatic procedure for cluster analysis using the baseline data
rng('default');  % set random seed for reproducibility
% load baseline data matrix
load('./results/dataBaseline.mat');
% step 1. data rescaling to 0 ~ 1 with min-max normalization, 22 variables
% are selected for rescaling. Gender is exclused from cluster analysis, so
% it's removed from the rescaling process.
[dataBaseline_rs,C,S] = normalize(table2array(dataBaseline(:,2:end)),'range'); 
% Step 2. use Calinski-Harabasz (C-H) psuedo F as a criteria to dermine the
% optimal number of clusters in the range of 2 to 10 clusters
    % Note: The Calinski-Harabasz index also known as the Variance Ratio Criterion, 
    % is the ratio of the sum of between-clusters dispersion and of inter-cluster dispersion for all clusters,
    % the higher the score , the better the performances.
eva = evalclusters(dataBaseline_rs,'kmeans','CalinskiHarabasz','KList',[2:10]); 
% find optimal cluster number K
K = eva.OptimalK;
% Plot the curve of C-H value against cluster number
figure; plot(eva); xlim([1 11]); pbaspect([2 1 1]); xlabel('Cluster number'); ylabel('C-H index'); title(strcat('Optimal cluster number is: ',num2str(K)));

% STEP 3. perform k-means clustering based on the optimal cluster number K,
% with 100 times of random initialization, and 100 iterations of cluster analysis
    % Note: The cluster analysis was repeated 100 times with random initialization to avoid final cluster membership derived from a local optimum. 
    % The maximal number of iterations was set at 100 to ensure ample iterations to arrive at an optimal clustering solution within each cluster analysis. 
[idx,cent,sumdist] = kmeans(dataBaseline_rs,K,'Distance','sqeuclidean', 'MaxIter',100,'Display','final','Replicates',100);                     

% save clustering results.
% The file contains cluster membership assignments. It also contains parameters needed as a decision rule 
% "idx" contains cluster memberships, 1.mild motor-non-motor subtype (MMNS)
%                                     2.severe motor-non-motor subtype (SMNS). 
% "cent" contains centroids for all 22 variables used in cluster analsys
% "sumdist" contains the sum of distance between samples and the centroid for each cluster
% "C" and "S" are the center and scaling factors for data rescaling based on
%             min-max normalization. The equation is data_rescaled = (data-C)/S
save('./results/clusteringResults','idx','cent','sumdist','C', 'S');
%% Part IV: Characterizing baseline variables based on PD subtype assignment
% add effect size toolbox to path, for calculating standardized mean difference (SMD)
addpath(genpath('./EffectSizeToolbox_v1.5'));
% creating a table with mean and std values for each variables for different subtypes
charactTable = array2table([mean(dataBaseline{idx==1,2:end}); std(dataBaseline{idx==1,2:end});
                            mean(dataBaseline{idx==2,2:end});std(dataBaseline{idx==2,2:end});],...
                            'VariableNames',dataBaseline.Properties.VariableNames(2:end));
charactTable.Properties.RowNames = {'MMNS Mean','MMNS SD','SMNS Mean','SMNS SD'};

% calculate standardized mean difference (SMD) between MMNS and SMNS using Hedges’ g stats.
pvals = []; smd = [];
for i = 2:size(dataBaseline,2) % first column is gender, which is not included in the statistical test
    pvals = [pvals ranksum(dataBaseline{idx==1,i},dataBaseline{idx==2,i})];
    smd = [smd; mes(dataBaseline{idx==1,i},dataBaseline{idx==2,i},'hedgesg','nBoot',3000)]; 
end
% adding SMD results into the characterization table
charactTable{5,:} = pvals; charactTable.Properties.RowNames{5} = 'P values';
charactTable{6,:} = [smd.hedgesg]; charactTable.Properties.RowNames{6} = 'SMD';
hedgesgCi=[smd.hedgesgCi];charactTable{7,:} = hedgesgCi(1,:);charactTable.Properties.RowNames{7} = 'loCI';
charactTable{8,:} = hedgesgCi(2,:); charactTable.Properties.RowNames{8} = 'hiCI';
% reorder the table
charactTable =rows2vars(charactTable);
charactTable_ro =charactTable([1 2 19:22 3:9 18 10:17],[1 2 4 3 5 6 7 8 9]);
% print out the table 
disp(charactTable_ro);

% visualize the SMD results for all 22 variables
figure;
for i = 1:size(charactTable_ro,1)
    plot([charactTable_ro{i,8},charactTable_ro{i,9}],[100-2*i,100-2*i],'k-');
    hold on;
    plot(charactTable_ro{i,7},100-2*i,'kd','MarkerSize',5);
end
hold off; xlim([-3 2]);
ylabel('22 variables');xlabel('SMD (MMNS - SMNS)');

% load baseline imaging variables and CSF variables
Data_bl = readtable('PPMI_all/PPMI_Baseline_Data_02Jul2018.csv');
imagingVarData = (Data_bl(ismember(Data_bl.PATNO,Baseline_allVar_rm.PATNO),85:112));
CSFVarData = (Data_bl(ismember(Data_bl.PATNO,Baseline_allVar_rm.PATNO),113:122));

% statistical test of imaging variables between PD subtypes
stats_imaging = [];
for i = 1:size(imagingVarData,2)
    group1 = imagingVarData{idx==1&~isnan(imagingVarData{:,i}),i};
    group2 = imagingVarData{idx==2&~isnan(imagingVarData{:,i}),i};
    stats_imaging(i,:) = [mean(group1) mean(group2) ranksum(group1,group2) std(group1) std(group2)];
end
stats_imaging_tb = array2table(stats_imaging,'VariableNames',{'MMNS mean', 'SMNS mean','p val','MMNS SD', 'SMNS SD'},'RowNames',imagingVarData.Properties.VariableNames);
disp(stats_imaging_tb);

% statistical test of CSF variables between PD subtypes
stats_CSF = [];
for i = 1:size(CSFVarData,2)
    group1 = CSFVarData{idx==1&~isnan(CSFVarData{:,i}),i};
    group2 = CSFVarData{idx==2&~isnan(CSFVarData{:,i}),i};
    stats_CSF(i,:) = [mean(group1) mean(group2) ranksum(group1,group2) std(group1) std(group2)];
end
stats_CSF_tb = array2table(stats_CSF,'VariableNames',{'MMNS mean', 'SMNS mean','p val','MMNS SD', 'SMNS SD'},'RowNames',CSFVarData.Properties.VariableNames);
disp(stats_CSF_tb);

% correlation between signficant CSF factors with MOCA and UPDRS III
[rho_moca,pval_moca] = corr(CSFVarData.ptau_asyn(~isnan(CSFVarData.ptau_asyn)),Baseline_allVar_rm.moca(~isnan(CSFVarData.ptau_asyn)));
[rho_UPDRS,pval_UPDRS] = corr(CSFVarData.ptau_asyn(~isnan(CSFVarData.ptau_asyn)),Baseline_allVar_rm.updrs3_score(~isnan(CSFVarData.ptau_asyn)));
fprintf('Correlation coefficent between MOCA and ptau_asyn ratio is %f, p = %f\n',rho_moca,pval_moca);
fprintf('Correlation coefficent between UPDRS III and ptau_asyn ratio is %f, p = %f\n',rho_UPDRS,pval_UPDRS);

% Calculate motor composite score, and non-motor composite score
motor_composite_score = mean(dataBaseline_rs(:,[3 20:23]),2);
nonmotor_composite_score = dataBaseline_rs(:,4:19);
% MOCA and UPSIT have oposite severity as other non-motor factors, and
% adjusted by substraction from 1
nonmotor_composite_score(:,1:2) = 1-nonmotor_composite_score(:,1:2);
nonmotor_composite_score = mean(nonmotor_composite_score,2);

% scatter plot of cluster membership along motor and non-motor axis.
figure;
gscatter(motor_composite_score,nonmotor_composite_score,idx,'rb',[],[20,20]);
hold on;
plot(mean(motor_composite_score(idx ==1)),mean(nonmotor_composite_score(idx ==1)),'k+','MarkerSize',15,'LineWidth',3);
plot(mean(motor_composite_score(idx ==2)),mean(nonmotor_composite_score(idx ==2)),'k+','MarkerSize',15,'LineWidth',3);
hold off;
legend({'MMNS','SMNS','Centroid'});
%% Part V: longitudinal evaluation
% load progression data, [ patno, H&Y year0-3, UPDRS-III year0-3, and MOCA year0-3]
load('./results/dataProgression.mat');
% load baseline data to retrieve patient ID with cluster memberships
load('./results/dataBaseline.mat','Baseline_allVar_rm');
% retrieve cluster memberships
load('./results/clusteringResults.mat','idx');

% remove pats with progression but not entered into cluster analysis due to missing key variables
dataProgression_rm=dataProgression(ismember(dataProgression(:,1),Baseline_allVar_rm.PATNO),:);

% adding cluster memberships to the last column
for i = 1:size(dataProgression_rm,1)
    dataProgression_rm(i,14)=idx(Baseline_allVar_rm.PATNO==dataProgression_rm(i,1));
end

% calcuate progression rate for early progression Y1-Y0, secondary Y3-Y1,
% and long-term progression Y3-Y0 for H&Y and MoCA
ProgRate = dataProgression_rm(:,[3 5 5 11 13 13])-dataProgression_rm(:,[2 3 2 10 11 10]);
Prog_mean = [nanmean(ProgRate(dataProgression_rm(:,14)==1,:))' nanmean(ProgRate(dataProgression_rm(:,14)==2,:))'];
Prog_SD = [nanstd(ProgRate(dataProgression_rm(:,14)==1,:))' nanstd(ProgRate(dataProgression_rm(:,14)==2,:))'];

ProgressionRateResults = array2table([Prog_mean Prog_SD],'VariableNames',{'MMNS mean', 'NMNS mean','MMNS SD', 'NMNS SD'},...
    'RowNames',{'Early Progression HY','Secondary Progression HY','Long-term Progression HY',...
    'Early Progression MOCA','Secondary Progression MOCA','Long-term Progression MOCA'});
disp(ProgressionRateResults);

% getting average H&Y and MoCA values for all 4 time points
MMNS_HY_mean_y0_3 = nanmean(dataProgression_rm(dataProgression_rm(:,14)==1,2:5),1);
MMNS_MOCA_mean_y0_3 = nanmean(dataProgression_rm(dataProgression_rm(:,14)==1,10:13),1);

SMNS_HY_mean_y0_3 = nanmean(dataProgression_rm(dataProgression_rm(:,14)==2,2:5),1);
SMNS_MOCA_mean_y0_3 = nanmean(dataProgression_rm(dataProgression_rm(:,14)==2,10:13),1);

% time plot of progression of both H&Y and MoCA
figure;
yyaxis left;
plot(0:3,MMNS_HY_mean_y0_3,'o-','MarkerSize',10);
hold on;
plot(0:3,SMNS_HY_mean_y0_3,'o--','MarkerSize',10);
hold off;
xlim([-0.5 3.5]);
ylabel('H&Y');
yyaxis right;
plot(0:3,MMNS_MOCA_mean_y0_3,'d-','MarkerSize',10);
hold on;
plot(0:3,SMNS_MOCA_mean_y0_3,'d--','MarkerSize',10);
hold off;
ylabel('MoCA');
legend({'MMNS H&Y','SMNS H&Y','MMNS MoCA','SMNS MoCA'});

% getting all H&Y and MoCA values for all 4 time points
MMNS_HY_y0_3 = dataProgression_rm(dataProgression_rm(:,14)==1,2:5);
MMNS_MOCA_y0_3 = dataProgression_rm(dataProgression_rm(:,14)==1,10:13);
SMNS_HY_y0_3 = dataProgression_rm(dataProgression_rm(:,14)==2,2:5);
SMNS_MOCA_y0_3 = dataProgression_rm(dataProgression_rm(:,14)==2,10:13);

% perform two-way ANOVA test to evaluate the effect PD subgroups and time points on MoCA and H&Y progression
% Two-way ANOVA for H&Y progression
[p_HY,tbl_HY,stats_HY] = anovan([MMNS_HY_y0_3(:);SMNS_HY_y0_3(:)],...
    {[ones(length(MMNS_HY_y0_3(:)),1)*1;ones(length(SMNS_HY_y0_3(:)),1)*2;],...
    [ones(size(MMNS_HY_y0_3,1),1)*1;ones(size(MMNS_HY_y0_3,1),1)*2;ones(size(MMNS_HY_y0_3,1),1)*3;ones(size(MMNS_HY_y0_3,1),1)*4;...
    ones(size(SMNS_HY_y0_3,1),1)*1;ones(size(SMNS_HY_y0_3,1),1)*2;ones(size(SMNS_HY_y0_3,1),1)*3;ones(size(SMNS_HY_y0_3,1),1)*4;]}...
    ,'model',2,'varnames',{'PD subtype','Time points'});
% Two-way ANOVA for MOCA progression
[p_MOCA,tbl_MOCA,stats_MOCA] = anovan([MMNS_MOCA_y0_3(:);SMNS_MOCA_y0_3(:)],...
    {[ones(length(MMNS_MOCA_y0_3(:)),1)*1;ones(length(SMNS_MOCA_y0_3(:)),1)*2;],...
    [ones(size(MMNS_MOCA_y0_3,1),1)*1;ones(size(MMNS_MOCA_y0_3,1),1)*2;ones(size(MMNS_MOCA_y0_3,1),1)*3;ones(size(MMNS_MOCA_y0_3,1),1)*4;...
    ones(size(SMNS_MOCA_y0_3,1),1)*1;ones(size(SMNS_MOCA_y0_3,1),1)*2;ones(size(SMNS_MOCA_y0_3,1),1)*3;ones(size(SMNS_MOCA_y0_3,1),1)*4;]}...
    ,'model',2,'varnames',{'PD subtype','Time points'});
