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
% Written by Jiwon Yeon, last edited Jan.23.2019.
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

% Accuracy
Observed = {round(mean(observed.answer2),3); []; []};
[h p ci stats] = ttest(observed.answer2, population.acc.answer2);
Population = {round(mean(population.acc.answer2),3); p; round(abs(stats.tstat),3)};
[h p ci stats] = ttest(observed.answer2, summary_random.acc.answer2);
Summary_Random = {round(mean(summary_random.acc.answer2),3); p; round(abs(stats.tstat),3)};
[h p ci stats] = ttest(observed.answer2, summary_strategic.acc.answer2);
Summary_Strategic = {round(mean(summary_strategic.acc.answer2),3); p; round(abs(stats.tstat),3)};

RowNames = {'Accuracy'; 't-val(vs. Observed)'; 'p-val(vs. Observed)'};
Prediction_2alternative = table(Observed,Population,Summary_Random,Summary_Strategic,...
    'RowNames', RowNames)


% AIC comparisons
for sub = 1:length(observed.answer2)
    AIC.population(sub) = population.resfit{sub}.AIC;    
    AIC.summary_random(sub) = summary_random.resfit{sub}.AIC;
    AIC.summary_strategic(sub) = summary_strategic.resfit{sub}.AIC;
end

% Compare AIC_Average
output = AICanalysis([mean(AIC.population) mean(AIC.summary_random)],'e');
output = [output; AICanalysis([mean(AIC.population) mean(AIC.summary_strategic)],'e')];    
Population = {output(1,1); output(2,1)};

output = AICanalysis([mean(AIC.summary_random) mean(AIC.summary_strategic)],'e');
Summary_Random = {[]; output(1,1)};

RowNames = {'Summary_Random vs.'; 'Summary_Strategic vs.'};
AIC_comparison_Average = table(Population, Summary_Random,'RowNames', RowNames)


% Compare AIC_Average
output = AICanalysis([sum(AIC.population) sum(AIC.summary_random)],'e');
output = [output; AICanalysis([sum(AIC.population) sum(AIC.summary_strategic)],'e')];    
Population = {output(1,1); output(2,1)};

output = AICanalysis([sum(AIC.summary_random) sum(AIC.summary_strategic)],'e');
Summary_Random = {[]; output(1,1)};

RowNames = {'Summary_Random vs.'; 'Summary_Strategic vs.'};
AIC_comparison_Total = table(Population, Summary_Random,'RowNames', RowNames)

