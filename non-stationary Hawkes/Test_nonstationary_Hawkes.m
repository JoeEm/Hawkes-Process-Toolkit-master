%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file for testing the non-stationary Hawkes processes 
% last modified data: 3/22/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%loading the data
load MiniBlogSeqsTwoAndHalfYears.mat
Data = MiniBlogSeqsTwoAndHalfYears(1);

alg1.LowRank = 0; % without low-rank regularizer
alg1.Sparse = 1; % with sparse regul  arizer
alg1.alphaS = 1;
alg1.GroupSparse = 1; % with group-sparse regularizer
alg1.alphaGS = 100;
alg1.outer = 8; %8
alg1.rho = 0.1; % the initial parameter for ADMM
alg1.inner = 5; %5
alg1.thres = 1e-5;
alg1.Tmax = [];
alg1.storeErr = 0;
alg1.storeLL = 0;
alg1.As = 10; %for Local Independence R.
alg1.Ap = 1000;%for Pairwise similarity R.

para = alg1;

FirstTimeFlag = 1;

%hyperparameter
Beta = 5000;
Alpha = 400;
ClusterNumbers = 2;
threshold = 500;




while(1)
    %Clustering Data
    ClusteredData = AssginPtsToCluster(Data, ClusterNumbers, para, FirstTimeFlag, Beta);
    FirstTimeFlag = 0;
 
    %Learning Hawkes Processes(MLE-SGL)for each Clustered-Data respectively
    for i =1:ClusterNumbers
        model(i) = Initialization_Basis(ClusteredData{i});
        model_MLE_S(i) = Learning_MLE_S_Basis(ClusteredData{i} , model(i), para);
    end
    
    
    %Comuputing the Loss
    Loss = TotalLoss( Alpha, Beta, HawkesLoglikelihood, ClusterNumbers, ClusteredData, Timeinterval)
    %it won't stop until stationarity of the parameters    
    if (FirstTimeFlag == 1)
        PrevLost = Loss;
    else
        if (abs(Loss - PrevLost)<threshold) %stationarity
            break;
        else
            PrevLost = Loss;
        end
    end
    %if the procedure doesn't stop, combine the CluterData again for convenience. 
    Totallength = 0;
    for g = 1:ClusterNumbers
       tmp = length(ClusteredData{g}.Time);
       Totallength = tmp + Totallength;
    end
    Data = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Feature', []);
        for i = 1:Totallength
                Data.Time(i) = Data.Time(i);
                Data.Mark(i) = Data.Mark(i);
                Data.Feature(i) = Data.Feature(i);
         end
end


%Visualization the Graph
for i = 1:ClusterNumbers
    figure 
         
    imagesc(model_MLE_S(i).A)
    title(sprintf('ClusterNum=%d ', i ));
end