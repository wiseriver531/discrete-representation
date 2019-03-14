# Discrete-representation

The analysis codes are written in MATLAB. Please refer to the below guidline:

@ System requirement
- The codes were run and tested on macOS (10.14.2) with MATLAB_R2016a

@ Installation guide
- Data and analysis codes for each experiment are saved in separate folders.
- To run analysis codes correctly, please don't move or rename any files or folders.
- Set the current directory in MATLAB to where the running code is located.

@ Instructions for use
- 'data' folder contains behavioral raw data ('subject_responses' folder) and model fitting outputs ('fitting results' folder). 
- 'behavioral_data_cleanup.m' file (dir: '/codes/analysis/organize behavioral responses') outputs 'dataForModeling.mat'. To run model fitting, 'dataForModeling.mat' file is needed.
- To check final results, run 'Run_code.m' file (dir: 'codes/analysis/print results/' folder in each experiment folder).
- To check instructions for reproduction, please read 'License' file. 

If you have any questions, write to Jiwon Yeon at j.yeon@gatech.edu or wiseriveer531@gmail.com
