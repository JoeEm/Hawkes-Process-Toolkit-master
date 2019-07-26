%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file for testing the KDD(TICC)-method Hawkes processes 
% last modified data: 7/12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear %clear workplace

%generating simulation data.
%'options' structure.
options.N = 500; % the number of sequences
options.Nmax = 500; % the maximum number of events per sequence
options.Tmax = 500; % the maximum size of time window
options.tstep = 0.1;
options.dt = 0.1; % the length of each time step
options.M = 250; % the number of steps in the time interval for computing sup-intensity
options.GenerationNum = 100; % the number of generations for branch processing
D = 4; % the dimension of Hawkes processes


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

para = alg1;

%hyperparameter

Beta = 1500;
IterationNum = 1;
FirstTimeLossFlag = 1;
% 
TotalIterationNum = 1000;
%modelStorage = cell(TotalIterationNum,1);
%Loss = zeros(TotalIterationNum,1);

%GenratingSimulationData

% MaskValue(1)=0.7;
% MaskValue(2)=0.7;
% amplifier(1) = 1;
% amplifier(2) = 4;
% ClusterNumbers = 4;
A1 = 0;
flag = 0;
ClusterNumbers = 4;
for j = 1:ClusterNumbers
   %amplifier = 1 + 9*(j-1);
    amplifier = 1;
    MaskValue = 0.7;
    if (j == 3)
       %A1 = paraDataGenrating{3}.A;
       A1 = paraDataGenrating{1}.A;
       flag = 1;
       group = 1;
       [Seqs1{j},paraDataGenrating{j}] = GenerateSingleGroupData(options,D,group,amplifier,MaskValue,A1, flag);
       Seqs2{j} = TailoringData(Seqs1{j});
    end
    if (j == 4 )
       A1 = paraDataGenrating{2}.A;
       flag = 1;
       group = 2;
       [Seqs1{j},paraDataGenrating{j}] = GenerateSingleGroupData(options,D,group,amplifier,MaskValue,A1, flag); 
       Seqs2{j} = TailoringData(Seqs1{j});
    end
    if (j == 1 || j ==2)
    %if (j == 1 || j ==2 || j == 3)
        [Seqs1{j},paraDataGenrating{j}] = GenerateSingleGroupData(options,D,j,amplifier,MaskValue,A1, flag);
        Seqs2{j} = TailoringData(Seqs1{j});
    end
    flag = 0;
end

%for more than 2 patterns 
for k = 1:(ClusterNumbers-1)
    Seqs2{k+1} = LinkingTwoSeqsToOne(Seqs2{k},Seqs2{k+1});
end
% Seqs = Seqs1{k+1};
% load Seqs3DifClu.mat
% load 3DifCluData29Iteration.mat
% load OriginalData_ABC4Dimen3000Seqs.mat
% load OriginalData_ABAB4Dimen1200Seqs.mat
% TmpData{1} = Data;
% Data = [];
% Seqs1 = TmpData{1};
% for k = 1:length(TmpData)
%    TmpData{k} = [];
% end

batchsize = 20; %parameter set mannually.
[NewSeqs,Data] = BatchSizingData(Seqs2{k+1} ,batchsize);
TmpData{1} =Data;

figure

ClusterNumbers = 2;
for IterationNum = 1:(TotalIterationNum-1)
    %Clustering Data
    fprintf('clustering Data\n');
    [NewSeqs,PointsLoss, location, FinalPath,Data] = FasterAssginPtsToCluster(Data, NewSeqs, ClusterNumbers, para, IterationNum, Beta,batchsize);
 
    %Comuputing the Loss
    fprintf('computing the Loss \n');
    SumLoss(IterationNum) = Loss(PointsLoss, Beta, FinalPath);
    TmpData{IterationNum + 1} = Data;
    %Visualization the Loss (and Graph)
    %Num = 1:1:TotalIterationNum-1;
    Num = 1:1:IterationNum;

    hold on
    plot(Num, SumLoss, 'bs-');
    axis tight
    axis square
    ylabel('Loss');
    xlabel('IterationNum');
   % hold off
    drawnow
end