%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file for testing the non-stationary Hawkes processes 
% last modified data: 3/22/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear %clear workplace
%loading the data
%load MiniBlogSeqsTwoAndHalfYears.mat

%preprecessing Data...(IPTVData...)
load IPTVData.mat %Seqs' length is 302 ,each 11 months long
%Seqs1 = Seqs(1:300);
Seqs1 = Seqs(1:30);
%NewSeqs = MiniBlogRealWorldCutting(MiniBlogSeqsTwoAndHalfYears);
NewSeqs2 = IPTV_OriginalRealData_Cutting(Seqs1);
NewSeqs1 = RearrangeTheMiniBlogSeqs(NewSeqs2);
NewSeqs = NewSeqs2(1:50);%extracting a years(50weeks aprox. to a year) data
NewSeqs = CombineTheOriginalData(NewSeqs,length(NewSeqs));
Data = NumberTheSeqs(NewSeqs); %add number label to each data


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

sample = 0;
 for k = 1:length(Data)
    sample = sample + length(Data(k).Time);
 end
sample = sample/length(Data);
alg1.N = sample;%number of sample

para = alg1;


%hyperparameter
Beta = 50;
Alpha = 400;
ClusterNumbers = 2;
threshold = 500;
Timeinterval = 1;

IterationNum = 1;
FirstTimeLossFlag = 1;
 

while(1)
    %Clustering Data
    fprintf('clustering Data\n');
    if  (IterationNum == 1)
        ClusteredData = AssginPtsToCluster(Data, ClusterNumbers, para, Beta,IterationNum);
    else
        ClusteredData = AssginPtsToCluster(NewData, ClusterNumbers, para, Beta,IterationNum ,model_MLE_S);
    end

    %dbstop in Test_nonstationary_Hawkes at 46;
    %Learning Hawkes Processes(MLE-SGL)for each Clustered-Data respectively
    fprintf('MLE-SGL learning\n');
    HawkesLoglikelihood = 0;
    for i =1:ClusterNumbers
        para.Tmax = ClusteredData{i}.Time(length(ClusteredData{i}.Time));
        %model(i) = Initialization_Basis(ClusteredData{i});
        %mannually set the parameters.
        model1(i).kernel = 'gauss';
        model1(i).w = pi*para.N/(para.Tmax);  %According to the Real-Data Experiment
        model1(i).landmark = (model1(i).w)*(0:ceil(para.Tmax/model1(i).w));
        D = 16; %dimension of IPTV data set.
        model1(i).A = rand(D, length(model1(i).landmark), D);
        model1(i).mu = rand(D, 1);
        %breking down the sequences
        OutSeqs1 = IPTV_RealData_Cutting(ClusteredData{i});
        OutSeqs = RearrangeTheMiniBlogSeqs(OutSeqs1);
        if (length(OutSeqs) == 1)
           model_MLE_S(i) = Learning_MLE_S_nonstationary(OutSeqs, model1(i), para);
        else
           model_MLE_S(i) = Learning_MLE_S_nonstationary(OutSeqs(1:length(OutSeqs)) , model1(i), para);
        end
        LogV = Loglike_Basis_NonStationary( OutSeqs(1:length(OutSeqs)), model_MLE_S(i), para );
        HawkesLoglikelihood = LogV + HawkesLoglikelihood;
    end
    
    
    %Comuputing the Loss
    fprintf('computing the Loss \n');
    Loss = TotalLoss( Alpha, Beta, HawkesLoglikelihood, ClusterNumbers, ClusteredData);
    %it won't stop until stationarity of the parameters
    if (FirstTimeLossFlag == 1)
        PrevLost = Loss;
    else
        if (abs(Loss - PrevLost)<threshold) %stationarity
            break;
        else
            PrevLost = Loss;
        end
    end
    FirstTimeLossFlag = 0;
    %if the procedure doesn't stop, combine the CluterData again for convenience. 
    NewData = CombineTheData(ClusteredData, ClusterNumbers);

    fprintf('abs(Loss - PrevLost) = %d, IterationNum =%d, \n',...
                abs(Loss - PrevLost), IterationNum); 
    IterationNum = IterationNum + 1; %counting
end


%Visualization the Graph
for i = 1:ClusterNumbers
    figure 
    A1 = ImpactFunc( model_MLE_S(i), para );

    imagesc(A1)
    title(sprintf('ClusterNum=%d ', i ));
end