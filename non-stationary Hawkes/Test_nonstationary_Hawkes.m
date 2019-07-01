%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file for testing the non-stationary Hawkes processes 
% last modified data: 3/22/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear %clear workplace

%loading the Real-World-Data
%load MiniBlogSeqsTwoAndHalfYears.mat
%preprecessing Data...(IPTVData...)
%load IPTVData.mat %Seqs' length is 302 ,each 11 months long
%Seqs1 = Seqs(1:300);
%Seqs1 = Seqs(1:30);
%NewSeqs = MiniBlogRealWorldCutting(MiniBlogSeqsTwoAndHalfYears);
%NewSeqs2 = IPTV_OriginalRealData_Cutting(Seqs1);
%NewSeqs1 = RearrangeTheMiniBlogSeqs(NewSeqs2);
%NewSeqs = NewSeqs2(1:50);%extracting a years(50weeks aprox. to a year) data
%NewSeqs = CombineTheOriginalData(NewSeqs,length(NewSeqs));
%Data = NumberTheSeqs(NewSeqs); %add number label to each data

%generating simulation data.
% 'options' structure.
options.N = 500; % the number of sequences
options.Nmax = 500; % the maximum number of events per sequence
options.Tmax = 200; % the maximum size of time window
options.tstep = 0.1;
options.dt = 0.1; % the length of each time step
options.M = 250; % the number of steps in the time interval for computing sup-intensity
options.GenerationNum = 100; % the number of generations for branch processing
D = 8; % the dimension of Hawkes processes


%BreakSeqs(InputSeqs,options)



alg1.LowRank = 0; % without low-rank regularizer
alg1.Sparse = 1; % with sparse regularizer
alg1.alphaS = 10; %MLE-S
alg1.GroupSparse = 1; % with group-sparse regularizer
alg1.alphaGS = 100; %MLE-SGL
alg1.outer = 5;%5
alg1.rho = 0.1; % the initial parameter for ADMM
alg1.inner = 8;%8
alg1.thres = 1e-5;
alg1.Tmax = [];
alg1.storeErr = 0;
alg1.storeLL = 0;
alg1.As = 10; %for Local Independence R.
alg1.Ap = 1000;%for Pairwise similarity R.

% sample = 0;
%  for k = 1:length(Data)
%     sample = sample + length(Data(k).Time);
%  end
% sample = sample/length(Data);
% alg1.N = sample;%number of sample

para = alg1;


%hyperparameter
ClusterNumbers = 2;

IterationNum = 1;
FirstTimeLossFlag = 1;

TotalIterationNum = 6;
%modelStorage = cell(TotalIterationNum,1);
%Loss = zeros(TotalIterationNum,1);

%GenratingSimulationData
for j = 1:ClusterNumbers
    [Seqs1{j},paraDataGenrating{j}] = GenerateSingleGroupData(options,D,j);
end
Seqs = LinkingTwoSeqsToOne(Seqs1{1},Seqs1{2});


for IterationNum = 1:TotalIterationNum
    %Clustering Data
    fprintf('clustering Data\n');
    [Outputmodel,ClusteredData] = AssginPtsToCluster(Seqs, ClusterNumbers, para, IterationNum, TotalIterationNum);
 
    %Comuputing the Loss
    fprintf('computing the Loss \n');
    Loss(IterationNum) = TotalLoss(ClusterNumbers, ClusteredData, Outputmodel, para);
    modelStorage{IterationNum} = Outputmodel;
end


%Visualization the Loss (and Graph)
Num = 1:1:TotalIterationNum;
figure

hold on
plot(Num, Loss, 'bs-');
axis tight
axis square
ylabel('Loss');
xlabel('Iteration');
hold off
