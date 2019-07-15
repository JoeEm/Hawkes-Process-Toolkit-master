%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file for testing the KDD(TICC)-method Hawkes processes 
% last modified data: 7/12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear %clear workplace

%generating simulation data.
%'options' structure.
options.N = 200; % the number of sequences
options.Nmax = 500; % the maximum number of events per sequence
options.Tmax = 500; % the maximum size of time window
options.tstep = 0.1;
options.dt = 0.1; % the length of each time step
options.M = 250; % the number of steps in the time interval for computing sup-intensity
options.GenerationNum = 1000; % the number of generations for branch processing
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
ClusterNumbers = 3;
Beta = 15000;
IterationNum = 1;
FirstTimeLossFlag = 1;

TotalIterationNum = 20;
%modelStorage = cell(TotalIterationNum,1);
%Loss = zeros(TotalIterationNum,1);

%GenratingSimulationData
for j = 1:ClusterNumbers
   %amplifier = 1 + 9*(j-1);
    amplifier = 1;
    [Seqs1{j},paraDataGenrating{j}] = GenerateSingleGroupData(options,D,j,amplifier);
end

%for more than 2 patterns 
for k = 1:(ClusterNumbers-1)
    Seqs1{k+1} = LinkingTwoSeqsToOne(Seqs1{k},Seqs1{k+1});
end
Seqs = Seqs1{k+1};
% load Seqs3DifClu.mat


batchsize = 25; %parameter set mannually.
NewSeqs = BatchSizingData(Seqs,batchsize);

figure

for IterationNum = 1:(TotalIterationNum-1)
    %Clustering Data
    fprintf('clustering Data\n');
    [NewSeqs,PointsLoss, location, FinalPath] = FasterAssginPtsToCluster(Seqs, NewSeqs, ClusterNumbers, para, IterationNum, Beta);
 
    %Comuputing the Loss
    fprintf('computing the Loss \n');
    SumLoss = Loss(PointsLoss, Beta, FinalPath);
    
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

%Visualization the Loss (and Graph)
%Num = 1:1:TotalIterationNum-1;
Num = 1:1:9;
figure

hold on
plot(Num, Loss, 'bs-');
axis tight
axis square
ylabel('Loss');
xlabel('Iteration');
hold off

figure

hold on
plot(Num, b, 'bs-');
axis tight
axis square
ylabel('Loglike of Part-2');
xlabel('Iteration');
hold off
