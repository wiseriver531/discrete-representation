# Discrete-representation

The analysis codes are written in MATLAB. Please refer to the below guidline:

@ System requirement
- The codes were run and tested on macOS (10.14.2) with MATLAB_R2016a

@ Installation guide
- Data and analysis codes for each experiment are saved in separate folders.
- To run analysis codes correctly, please don't move or rename any files or folders.
- Set the current directory in MATLAB to where the running code is located.

@ Instructions for use
- Three experiments data are saved in separate folders.
- Under each experiment folder, 'codes' folder contains m-files for running model fitting and 'data' folder contains behavioral raw data ('subject_responses' folder) and model fitting outputs ('fitting results' folder). 
- 'behavioral_data_cleanup.m' file (dir: '/codes/analysis/organize behavioral responses/') outputs 'dataForModeling.mat'. To run model fitting, 'dataForModeling.mat' file is needed.
- To check final outputs from already saved fitting results, run 'Run_code.m' file (dir: '/codes/analysis/print results/'.
- To check instructions for reproduction, please read 'License' file. 

@ Demo
- To check the model fitting process from the beginning, first make sure to have 'dataForModeling.mat' file under '/data/subject_responses/' directory
- 'fitting_extended' folder is for the model that has more number of parameters (Exp1: 12; Exp2&3: 30).
- 'fitting_simple' folder is for the model that has reduced number of parameters (Exp1: 8; Exp2&3: 12).
- Run codes in order: 'step1_initial_parameter_fitting.m', 'step2_additional_parameter_fitting.m', and 'step3_test_fitting.m' (dir: '/codes/analysis/fitting_*/').
- Step1 code will give initial parameter fitting result.
- Step2 code will generate more parameter fitting result using different initial parameter set.
- The entire procedure will take a few hours to run on a "normal" computer.

If you have any questions, write to Jiwon Yeon at j.yeon@gatech.edu or wiseriveer531@gmail.com
