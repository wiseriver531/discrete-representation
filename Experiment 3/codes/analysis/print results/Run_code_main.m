% ------------------------------------------------------------------------
% Data analysis code for the manuscript "The nature of the perceptual 
% representation for decision making".
%
% To run this code, locate your current directory to where the code is 
% saved. If you want to see the results of the main analysis (i.e., with 
% extended parameters) type 'extended' for the 'version' variable. If you 
% want to check the results of the simple version analysis (i.e., with
% fewer parameters) type 'simple' for the 'version' variable. 
%
% Written by Jiwon Yeon, last edited Sep.14.2019.
% ------------------------------------------------------------------------
clear all, clc
version = 'extended';   % 'extended' or 'simple'

%% load, organize, test, and print data
dataPath = fileparts(fileparts(fileparts(pwd)));
load([dataPath '/data/subject_responses/dataForModeling.mat'])
observed.answer1 = data.acc(:,1);
observed.answer2 = data.acc(:,2);

load([dataPath '/data/fitting results/' version '/population_' version '.mat'])
population.acc.answer1 = accuracy_cond1;
population.acc.answer2 = accuracy_cond2;
population.resfit = resfit;

load([dataPath '/data/fitting results/' version '/summary+random_' version '.mat'])
summary_random.acc.answer1 = accuracy_cond1;
summary_random.acc.answer2 = accuracy_cond2;
summary_random.resfit = resfit;

load([dataPath '/data/fitting results/' version '/summary+strategic_' version '.mat'])
summary_strategic.acc.answer1 = accuracy_cond1;
summary_strategic.acc.answer2 = accuracy_cond2;
summary_strategic.resfit = resfit;


%% Accuracy
% 1st choice
RowNames = {'Accuracy'};
Observed = {round(mean(observed.answer1),3)};
Population = {round(mean(population.acc.answer1),3)};
Summary_Random = {round(mean(summary_random.acc.answer1),3)};
Summary_Strategic = {round(mean(summary_strategic.acc.answer1),3)};
Prediction_1st_answer = table(Observed,Population,Summary_Random,Summary_Strategic,...
    'RowNames', RowNames)

% 2nd choice
RowNames = {'Accuracy'; 'p-val(vs. Observed)'; 't-val(vs. Observed)'};
Observed = {round(mean(observed.answer2),3); []; []};
[h p ci stats] = ttest(observed.answer2, population.acc.answer2);
Population = {round(mean(population.acc.answer2),3); p; abs(round((stats.tstat),3))};
[h p ci stats] = ttest(observed.answer2, summary_random.acc.answer2);
Summary_Random = {round(mean(summary_random.acc.answer2),3); p; abs(round((stats.tstat),3))};
[h p ci stats] = ttest(observed.answer2, summary_strategic.acc.answer2);
Summary_Strategic = {round(mean(summary_strategic.acc.answer2),3); p; abs(round((stats.tstat),3))};
Prediction_2nd_answer = table(Observed,Population,Summary_Random,Summary_Strategic,...
    'RowNames', RowNames)


%% AIC comparisons
for sub = 1:length(observed.answer2)
    AIC.population(sub) = population.resfit{sub}.AIC;    
    AIC.summary_random(sub) = summary_random.resfit{sub}.AIC;
    AIC.summary_strategic(sub) = summary_strategic.resfit{sub}.AIC;
end

output = AICanalysis([mean(AIC.population) mean(AIC.summary_random)],'e');
Average_sum_random = [mean(AIC.population)-mean(AIC.summary_random);output(1,1)];
output = AICanalysis([sum(AIC.population) sum(AIC.summary_random)],'e');
Total_sum_random = [sum(AIC.population)-sum(AIC.summary_random); output(1,1)];

output = AICanalysis([mean(AIC.population) mean(AIC.summary_strategic)],'e');
Average_sum_strategic = [mean(AIC.population)-mean(AIC.summary_strategic);output(1,1)];
output = AICanalysis([sum(AIC.population) sum(AIC.summary_strategic)],'e');
Total_sum_strategic = [sum(AIC.population)-sum(AIC.summary_strategic); output(1,1)];

AIC_comparison_Population_vs_Summary = table(Average_sum_random, Total_sum_random, Average_sum_strategic, ...
    Total_sum_strategic, 'RowNames', {'Difference', 'Evidence ratio'})
