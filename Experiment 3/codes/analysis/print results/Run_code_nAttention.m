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
observed.answer1 = data.acc(:,1);
observed.answer2 = data.acc(:,2);

version = 'extended'; 
load([dataPath '/data/fitting results/' version '/summary+random_' version '.mat'])
summary_random.acc.answer1 = accuracy_cond1;
summary_random.acc.answer2 = accuracy_cond2;
summary_random.resfit = resfit;

load([dataPath '/data/fitting results/attention_extended/attention_2.mat'])
att2.acc.answer1 = accuracy_cond1;
att2.acc.answer2 = accuracy_cond2;
att2.resfit = resfit;

load([dataPath '/data/fitting results/attention_extended/attention_3.mat'])
att3.acc.answer1 = accuracy_cond1;
att3.acc.answer2 = accuracy_cond2;
att3.resfit = resfit;


%% Accuracy
% 1st choice
RowNames = {'Accuracy'; 'p-val(vs. Observed)'; 't-val(vs. Observed)'};
Observed = {round(mean(observed.answer1),3); []; []};
Summary_Random = {round(mean(summary_random.acc.answer1),3); []; []};
[h p ci stats] = ttest(observed.answer1, att2.acc.answer1);
Attention2 = {round(mean(att2.acc.answer1),3); p; abs(round((stats.tstat),3))};
[h p ci stats] = ttest(observed.answer1, att3.acc.answer1);
Attention3 = {round(mean(att3.acc.answer1),3); p; abs(round((stats.tstat),3))};
Prediction_1st_answer = table(Observed,Summary_Random,Attention2,Attention3,...
    'RowNames', RowNames)

% 2nd choice
RowNames = {'Accuracy'; 'p-val(vs. Observed)'; 't-val(vs. Observed)'};
Observed = {round(mean(observed.answer2),3); []; []};
Summary_Random = {round(mean(summary_random.acc.answer2),3); []; []};
[h p ci stats] = ttest(observed.answer2, att2.acc.answer2);
Attention2 = {round(mean(att2.acc.answer2),3); p; abs(round((stats.tstat),3))};
[h p ci stats] = ttest(observed.answer2, att3.acc.answer2);
Attention3 = {round(mean(att3.acc.answer2),3); p; abs(round((stats.tstat),3))};
Prediction_1st_answer = table(Observed,Summary_Random,Attention2,Attention3,...
    'RowNames', RowNames)

%% AIC comparisons
for sub = 1:length(observed.answer2)
    AIC.summary_random(sub) = summary_random.resfit{sub}.AIC;
    AIC.att2(sub) = att2.resfit{sub}.AIC;
    AIC.att3(sub) = att3.resfit{sub}.AIC;
end

output = AICanalysis([mean(AIC.att2) mean(AIC.summary_random)],'e');
Average_att2 = [mean(AIC.att2)-mean(AIC.summary_random); output(1,1)];

output = AICanalysis([mean(AIC.att3) mean(AIC.summary_random)],'e');
Average_att3 = [mean(AIC.att3)-mean(AIC.summary_random); output(1,1)];

output = AICanalysis([sum(AIC.att2) sum(AIC.summary_random)],'e');
Total_att2 = [sum(AIC.att2)-sum(AIC.summary_random); output(1,1)];

output = AICanalysis([sum(AIC.att3) sum(AIC.summary_random)],'e');
Total_att3 = [sum(AIC.att3)-sum(AIC.summary_random); output(1,1)];

AIC_comparison_SummaryRandom_vs_nAttention = table(Average_att2,Average_att3,Total_att2,Total_att3,...
    'RowNames', {'Difference', 'Evidence ratio'})
