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
observed.alternative6 = data.acc(:,1);
observed.alternative2 = data.acc(:,2);

load([dataPath '/data/fitting results/' version '/population_' version '.mat'])
population.acc.alternative6 = accuracy_cond1;
population.acc.alternative2 = accuracy_cond2;
population.resfit = resfit;

load([dataPath '/data/fitting results/' version '/summary_' version '.mat'])
summary.acc.alternative6 = accuracy_cond1;
summary.acc.alternative2 = accuracy_cond2;
summary.resfit = resfit;

%% Accuracy
% 6-alternative
RowNames = {'Accuracy'};
Observed = {round(mean(observed.alternative6),3)};
Population = {round(mean(population.acc.alternative6),3)};
Summary = {round(mean(summary.acc.alternative6),3)};
Prediction_6alternative = table(Observed,Population,Summary,...
    'RowNames', RowNames)

% 2-alternatvie
RowNames = {'Accuracy'; 'p-val(vs. Observed)'; 't-val(vs. Observed)'};
Observed = {round(mean(observed.alternative2),3); []; []};
[h p ci stats] = ttest(observed.alternative2, population.acc.alternative2);
Population = {round(mean(population.acc.alternative2),3); p; round(abs(stats.tstat),3)};
[h p ci stats] = ttest(observed.alternative2, summary.acc.alternative2);
Summary = {round(mean(summary.acc.alternative2),3); p; round(abs(stats.tstat),3)};
Prediction_2alternative = table(Observed,Population,Summary,...
    'RowNames', RowNames)


%% AIC comparisons
for sub = 1:length(observed.alternative2)
    AIC.population(sub) = population.resfit{sub}.AIC;    
    AIC.summary(sub) = summary.resfit{sub}.AIC; 
end

output = AICanalysis([mean(AIC.population) mean(AIC.summary)],'e');
Average = [mean(AIC.population)-mean(AIC.summary);output(1,1)];
output = AICanalysis([sum(AIC.population) sum(AIC.summary)],'e');
Total = [sum(AIC.population)-sum(AIC.summary); output(1,1)];
AIC_comparison_Population_vs_Summary = table(Average, Total, 'RowNames', {'Difference', 'Evidence ratio'})
