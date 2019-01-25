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
% Written by Jiwon Yeon, last edited Jan.22.2019.
% ------------------------------------------------------------------------
clear all, clc
version = 'extended';   % 'extended' or 'simple'


%% load, organize, test, and print data
folderName = ['fitting results/modeling_' version];
dataPath = [fileparts(fileparts(fileparts(pwd))) '/data/'];
load([dataPath '/subject_responses/dataForModeling'])
observed.alternative4 = data.accuracy(:,1);
observed.alternative2 = data.accuracy(:,2);

load([dataPath folderName '/population_' version '.mat'])
population.acc.alternative4 = accuracy_cond1;
population.acc.alternative2 = accuracy_cond2;
population.resfit = resfit;

load([dataPath folderName '/summary_' version '.mat'])
summary.acc.alternative4 = accuracy_cond1;
summary.acc.alternative2 = accuracy_cond2;
summary.resfit = resfit;

load([dataPath folderName '/twohighest_' version '.mat'])
twohighest.acc.alternative4 = accuracy_cond1;
twohighest.acc.alternative2 = accuracy_cond2;
twohighest.resfit = resfit;

% Accuracy
Observed = {round(mean(observed.alternative2),3); []; []};
[h p ci stats] = ttest(observed.alternative2, population.acc.alternative2);
Population = {round(mean(population.acc.alternative2),3); p; round(abs(stats.tstat),3)};
[h p ci stats] = ttest(observed.alternative2, summary.acc.alternative2);
Summary = {round(mean(summary.acc.alternative2),3); p; round(abs(stats.tstat),3)};
[h p ci stats] = ttest(observed.alternative2, twohighest.acc.alternative2);
TwoHighest = {round(mean(twohighest.acc.alternative2),3); p; round(abs(stats.tstat),3)};

RowNames = {'Accuracy'; 't-val(vs. Observed)'; 'p-val(vs. Observed)'};
Prediction_2alternative = table(Observed,Population,Summary,TwoHighest,...
    'RowNames', RowNames)


% AIC comparisons
for sub = 1:length(observed.alternative2)
    AIC.population(sub) = population.resfit{sub}.AIC;    
    AIC.summary(sub) = summary.resfit{sub}.AIC;
    AIC.twohighest(sub) = twohighest.resfit{sub}.AIC;
end

% Compare AIC_Average
output = AICanalysis([mean(AIC.population) mean(AIC.summary)],'e');
output = [output; AICanalysis([mean(AIC.population) mean(AIC.twohighest)],'e')];
Population = {output(1,1); output(2,1)};

output = AICanalysis([mean(AIC.twohighest) mean(AIC.summary)],'e');
TwoHighest = {output(1,1); []};

RowNames = {'Summary vs.'; 'TwoHighest vs.'};
AIC_comparison_Average = table(Population, TwoHighest, 'RowNames', RowNames)


% Compare AIC_Total
output = AICanalysis([sum(AIC.population) sum(AIC.summary)],'e');
output = [output; AICanalysis([sum(AIC.population) sum(AIC.twohighest)],'e')];
Population = {output(1,1); output(2,1)};

output = AICanalysis([sum(AIC.twohighest) sum(AIC.summary)],'e');
TwoHighest = {output(1,1); []};

RowNames = {'Summary vs.'; 'TwoHighest vs.'};
AIC_comparison_Total = table(Population, TwoHighest, 'RowNames', RowNames)
