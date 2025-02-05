%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file for testing the non-stationary Hawkes processes (Sythetic Data)
% last modified data: 3/22/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear %clear workplace

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
alg1.Sparse=1;
alg1.GroupSparse=1; 

para = alg1;

%generating SytheticData(containing Two pattarns)
Data1 = DataGeneratingNonStationary(para);%1-2000/2001-4000
Data = NumberTheSeqs(Data1);%add number label for each Seqs, and also intialize the 'Feature'

%hyperparameter
Beta = 5;
Alpha = 400;
ClusterNumbers = 2;
threshold = 500;

IterationNum = 1;
FirstTimeLossFlag = 1;
 

while(1)
    %Clustering Data
    fprintf('clustering Data\n');
    if  (IterationNum == 1)
        [OutputModel,ClusteredData] = AssginPtsToCluster(Data, ClusterNumbers, para, Beta,IterationNum);
    else
        [OutputModel,ClusteredData] = AssginPtsToCluster(NewData, ClusterNumbers, para, Beta,IterationNum ,OutputModel);
    end

    %dbstop in Test_nonstationary_Hawkes at 46;
    %Learning Hawkes Processes(MLE-SGL)for each Clustered-Data respectively
    fprintf('MLE-SGL learning\n');
    HawkesLoglikelihood = 0;
    for k =1:ClusterNumbers
        LogV = Loglike_Basis_NonStationary(ClusteredData{k},  OutputModel(k), para);
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
    %NewData = CombineTheData(ClusteredData, ClusterNumbers);
    NewData = ClusteredData;
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