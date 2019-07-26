%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Setup THAP: the toolkit for learning and analysis of Hawkes processes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

addpath('BasicFunc');
addpath('Data');
addpath('Simulation');
addpath('Learning');
addpath('Analysis');
addpath('Visualization');
addpath('MiniBlogData');
addpath('non-stationary Hawkes');
addpath('MCMCLearning');
addpath('NonStationarySeqsTestData');
addpath('resultAofMCMC');
addpath('KDDHawkes');
addpath('experiment720');

load('MiniBlogSeqsTwoAndHalfYears.mat');
Seqs1 = MiniBlogSeqsTwoAndHalfYears;
for k = 1:length(MiniBlogSeqsTwoAndHalfYears(2).Time)
   MiniBlogSeqsTwoAndHalfYears(2).Time(k) = MiniBlogSeqsTwoAndHalfYears(2).Time(k) + MiniBlogSeqsTwoAndHalfYears(1).Time(end);
end