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

load([dataPath '/data/fitting results/' version '/twohighest_' version '.mat'])
twohighest.acc.alternative6 = accuracy_cond1;
twohighest.acc.alternative2 = accuracy_cond2;
twohighest.resfit = resfit;

load([dataPath '/data/fitting results/' version '/threehighest_' version '.mat'])
threehighest.acc.alternative6 = accuracy_cond1;
threehighest.acc.alternative2 = accuracy_cond2;
threehighest.resfit = resfit;

% Accuracy
Observed = {round(mean(observed.alternative2),3); []; []};
[h p ci stats] = ttest(observed.alternative2, population.acc.alternative2);
Population = {round(mean(population.acc.alternative2),3); p; round(abs(stats.tstat),3)};
[h p ci stats] = ttest(observed.alternative2, summary.acc.alternative2);
Summary = {round(mean(summary.acc.alternative2),3); p; round(abs(stats.tstat),3)};
[h p ci stats] = ttest(observed.alternative2, twohighest.acc.alternative2);
TwoHighest = {round(mean(twohighest.acc.alternative2),3); p; round(abs(stats.tstat),3)};
[h p ci stats] = ttest(observed.alternative2, threehighest.acc.alternative2);
ThreeHighest = {round(mean(threehighest.acc.alternative2),3); p; round(abs(stats.tstat),3)};

RowNames = {'Accuracy'; 't-val(vs. Observed)'; 'p-val(vs. Observed)'};
Prediction_2alternative = table(Observed,Population,Summary,TwoHighest,ThreeHighest,...
    'RowNames', RowNames)


% AIC comparisons
for sub = 1:length(observed.alternative2)
    AIC.population(sub) = population.resfit{sub}.AIC;    
    AIC.summary(sub) = summary.resfit{sub}.AIC;
    AIC.twohighest(sub) = twohighest.resfit{sub}.AIC;
    AIC.threehighest(sub) = threehighest.resfit{sub}.AIC;
end

% Compare AIC_Average
output = AICanalysis([mean(AIC.population) mean(AIC.summary)],'e');
output = [output; AICanalysis([mean(AIC.population) mean(AIC.threehighest)],'e')];
output = [output; AICanalysis([mean(AIC.population) mean(AIC.twohighest)],'e')];
Population = {output(1,1); output(2,1); output(3,1)};

output = AICanalysis([mean(AIC.threehighest) mean(AIC.summary)],'e');
ThreeHighest = {output(1,1); []; []};

output = AICanalysis([mean(AIC.twohighest) mean(AIC.summary)],'e');
TwoHighest = {output(1,1); []; []};

RowNames = {'Summary vs.'; 'ThreeHighest vs'; 'TwoHighest vs.'};
AIC_comparison_Average = table(Population, ThreeHighest, TwoHighest, 'RowNames', RowNames)


% Compare AIC_Total
output = AICanalysis([sum(AIC.population) sum(AIC.summary)],'e');
output = [output; AICanalysis([sum(AIC.population) sum(AIC.threehighest)],'e')];
output = [output; AICanalysis([sum(AIC.population) sum(AIC.twohighest)],'e')];
Population = {output(1,1); output(2,1); output(3,1)};

output = AICanalysis([sum(AIC.threehighest) sum(AIC.summary)],'e');
ThreeHighest = {output(1,1); []; []};

output = AICanalysis([sum(AIC.twohighest) sum(AIC.summary)],'e');
TwoHighest = {output(1,1); []; []};

RowNames = {'Summary vs.'; 'ThreeHighest vs'; 'TwoHighest vs.'};
AIC_comparison_Total = table(Population, ThreeHighest, TwoHighest, 'RowNames', RowNames)
