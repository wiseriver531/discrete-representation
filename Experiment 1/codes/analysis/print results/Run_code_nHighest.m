% ------------------------------------------------------------------------
% Data analysis code for the manuscript "The nature of the perceptual 
% representation for decision making".
%
% To run this code, locate your current directory to where the code is 
% saved. The code is for result comparsions between summary and n-Highest 
% models. 
%
% Written by Jiwon Yeon, last edited Sep.14.2019.
% ------------------------------------------------------------------------
clear all, clc

%% load, organize, test, and print data
dataPath = [fileparts(fileparts(fileparts(pwd))) '/data/'];
load([dataPath '/subject_responses/dataForModeling'])
observed.alternative4 = data.accuracy(:,1);
observed.alternative2 = data.accuracy(:,2);

load([dataPath  '/fitting results/extended/summary_extended.mat'])
summary.acc.alternative4 = accuracy_cond1;
summary.acc.alternative2 = accuracy_cond2;
summary.resfit = resfit;

load([dataPath '/fitting results/extended/twohighest_extended.mat'])
high2.acc.alternative4 = accuracy_cond1;
high2.acc.alternative2 = accuracy_cond2;
high2.resfit = resfit;


%% Accuracy
% 4-alternative
RowNames = {'Accuracy'};
Observed = {round(mean(observed.alternative4),3)};
Summary = {round(mean(summary.acc.alternative4),3);};
Highest2 = {round(mean(high2.acc.alternative4),3)};
Prediction_6alternative = table(Observed,Summary,Highest2,...
    'RowNames', RowNames)

% 2-alternative
RowNames = {'Accuracy'; 'p-val(vs. Observed)'; 't-val(vs. Observed)'};
Observed = {round(mean(observed.alternative2),3); []; []};
[h p ci stats] = ttest(observed.alternative2, summary.acc.alternative2);
Summary = {round(mean(summary.acc.alternative2),3); p; round(abs(stats.tstat),3)};
[h p ci stats] = ttest(observed.alternative2, high2.acc.alternative2);
Highest2 = {round(mean(high2.acc.alternative2),3); p; round(abs(stats.tstat),3)};
Prediction_2alternative = table(Observed,Summary,Highest2,...
    'RowNames', RowNames)

%% AIC comparisons
for sub = 1:length(observed.alternative2)
    AIC.summary(sub) = summary.resfit{sub}.AIC;
    AIC.high2(sub) = high2.resfit{sub}.AIC;
end

output = AICanalysis([mean(AIC.high2) mean(AIC.summary)],'e');
Average_high2 = [mean(AIC.high2)-mean(AIC.summary); output(1,1)];

output = AICanalysis([sum(AIC.high2) sum(AIC.summary)],'e');
Total_high2 = [sum(AIC.high2)-sum(AIC.summary); output(1,1)];

AIC_comparison_Summary_vs_nHighest = table(Average_high2,Total_high2,...
    'RowNames', {'Difference', 'Evidence ratio'})
