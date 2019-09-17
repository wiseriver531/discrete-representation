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
dataPath = fileparts(fileparts(fileparts(pwd)));
load([dataPath '/data/subject_responses/dataForModeling.mat'])
observed.alternative6 = data.acc(:,1);
observed.alternative2 = data.acc(:,2);

version = 'extended';  
load([dataPath '/data/fitting results/' version '/summary_' version '.mat'])
summary.acc.alternative6 = accuracy_cond1;
summary.acc.alternative2 = accuracy_cond2;
summary.resfit = resfit;

load([dataPath '/data/fitting results/' version '/twohighest_' version '.mat'])
high2.acc.alternative6 = accuracy_cond1;
high2.acc.alternative2 = accuracy_cond2;
high2.resfit = resfit;

load([dataPath '/data/fitting results/' version '/threehighest_' version '.mat'])
high3.acc.alternative6 = accuracy_cond1;
high3.acc.alternative2 = accuracy_cond2;
high3.resfit = resfit;

%% Accuracy
% 6-alternative
RowNames = {'Accuracy'};
Observed = {round(mean(observed.alternative6),3)};
Summary = {round(mean(summary.acc.alternative6),3)};
Highest2 = {round(mean(high2.acc.alternative6),3)};
Highest3 = {round(mean(high3.acc.alternative6),3)};
Prediction_6alternative = table(Observed,Summary,Highest2,Highest3,...
    'RowNames', RowNames)

% 2-alternatvie
RowNames = {'Accuracy'; 'p-val(vs. Observed)'; 't-val(vs. Observed)'};
Observed = {round(mean(observed.alternative2),3); []; []};
[h p ci stats] = ttest(observed.alternative2, summary.acc.alternative2);
Summary = {round(mean(summary.acc.alternative2),3); p; round(abs(stats.tstat),3)};
[h p ci stats] = ttest(observed.alternative2, high2.acc.alternative2);
Highest2 = {round(mean(high2.acc.alternative2),3); p; round(abs(stats.tstat),3)};
[h p ci stats] = ttest(observed.alternative2, high3.acc.alternative2);
Highest3 = {round(mean(high3.acc.alternative2),3); p; round(abs(stats.tstat),3)};
Prediction_6alternative = table(Observed,Summary,Highest2,Highest3,...
    'RowNames', RowNames)


%% AIC comparisons
for sub = 1:length(observed.alternative2)
    AIC.summary(sub) = summary.resfit{sub}.AIC;
    AIC.high2(sub) = high2.resfit{sub}.AIC;
    AIC.high3(sub) = high3.resfit{sub}.AIC;
end

output = AICanalysis([mean(AIC.high2) mean(AIC.summary)],'e');
Average_high2 = [mean(AIC.high2)-mean(AIC.summary); output(1,1)];

output = AICanalysis([mean(AIC.high3) mean(AIC.summary)],'e');
Average_high3 = [mean(AIC.high3)-mean(AIC.summary); output(1,1)];

output = AICanalysis([sum(AIC.high2) sum(AIC.summary)],'e');
Total_high2 = [sum(AIC.high2)-sum(AIC.summary); output(1,1)];

output = AICanalysis([sum(AIC.high3) sum(AIC.summary)],'e');
Total_high3 = [sum(AIC.high3)-sum(AIC.summary); output(1,1)];

AIC_comparison_Summary_vs_nHighest = table(Average_high2,Average_high3,Total_high2,Total_high3,...
    'RowNames', {'Difference', 'Evidence ratio'})
