% ------------------------------------------------------------------------
% Data analysis code for the manuscript "The nature of the perceptual 
% representation for decision making".
%
% To run this code, locate your current directory to where the code is 
% saved. The code is for result comparsions between summary and n-Attention 
% models. 
%
% Written by Jiwon Yeon, last edited Sep.14.2019.
% ------------------------------------------------------------------------
clear all, clc

%% load, organize, test, and print data
dataPath = fileparts(fileparts(fileparts(pwd)));
load([dataPath '/data/subject_responses/dataForModeling.mat'])
observed.alternative6 = data.acc(:,1);
observed.alternative2 = data.acc(:,2);

version = 'extended';  
load([dataPath '/data/fitting results/' version '/summary_' version '.mat'])
summary.acc.alternative6 = accuracy_cond1;
summary.acc.alternative2 = accuracy_cond2;
summary.resfit = resfit;

load([dataPath '/data/fitting results/attention_extended/attention_2.mat'])
att2.acc.alternative6 = accuracy_cond1;
att2.acc.alternative2 = accuracy_cond2;
att2.resfit = resfit;

load([dataPath '/data/fitting results/attention_extended/attention_3.mat'])
att3.acc.alternative6 = accuracy_cond1;
att3.acc.alternative2 = accuracy_cond2;
att3.resfit = resfit;


%% Accuracy
% 6-alternative
RowNames = {'Accuracy'; 'p-val(vs. Observed)'; 't-val(vs. Observed)'};
Observed = {round(mean(observed.alternative6),3); []; []};
[h p ci stats] = ttest(observed.alternative6, summary.acc.alternative6);
Summary = {round(mean(summary.acc.alternative6),3); p; round(abs(stats.tstat),3)};
[h p ci stats] = ttest(observed.alternative6, att2.acc.alternative6);
Attention2 = {round(mean(att2.acc.alternative6),3); p; round(abs(stats.tstat),3)};
[h p ci stats] = ttest(observed.alternative6, att3.acc.alternative6);
Attention3 = {round(mean(att3.acc.alternative6),3); p; round(abs(stats.tstat),3)};
Prediction_6alternative = table(Observed,Summary,Attention2,Attention3,...
    'RowNames', RowNames)

% 2-alternatvie
Observed = {round(mean(observed.alternative2),3); []; []};
[h p ci stats] = ttest(observed.alternative2, summary.acc.alternative2);
Summary = {round(mean(summary.acc.alternative2),3); p; round(abs(stats.tstat),3)};
[h p ci stats] = ttest(observed.alternative2, att2.acc.alternative2);
Attention2 = {round(mean(att2.acc.alternative2),3); p; round(abs(stats.tstat),3)};
[h p ci stats] = ttest(observed.alternative2, att3.acc.alternative2);
Attention3 = {round(mean(att3.acc.alternative2),3); p; round(abs(stats.tstat),3)};
Prediction_6alternative = table(Observed,Summary,Attention2,Attention3,...
    'RowNames', RowNames)

%% AIC comparisons
for sub = 1:length(observed.alternative2)
    AIC.summary(sub) = summary.resfit{sub}.AIC;
    AIC.att2(sub) = att2.resfit{sub}.AIC;
    AIC.att3(sub) = att3.resfit{sub}.AIC;
end

output = AICanalysis([mean(AIC.att2) mean(AIC.summary)],'e');
Average_att2 = [mean(AIC.att2)-mean(AIC.summary); output(1,1)];

output = AICanalysis([mean(AIC.att3) mean(AIC.summary)],'e');
Average_att3 = [mean(AIC.att3)-mean(AIC.summary); output(1,1)];

output = AICanalysis([sum(AIC.att2) sum(AIC.summary)],'e');
Total_att2 = [sum(AIC.att2)-sum(AIC.summary); output(1,1)];

output = AICanalysis([sum(AIC.att3) sum(AIC.summary)],'e');
Total_att3 = [sum(AIC.att3)-sum(AIC.summary); output(1,1)];

AIC_comparison_Summary_vs_nAttention = table(Average_att2,Average_att3,Total_att2,Total_att3,...
    'RowNames', {'Difference', 'Evidence ratio'})
