%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AssginPtsToCluster() : The function receives the orignial Data , 
% then Assgins the points to corresponding cluster via Algorithnm 1,
% and finally return the ClusteredData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Outputmodel,ClusteredData] = AssginPtsToCluster(Data, ClusterNumbers, para, IterationNum, TotalIterationNum)
%Initialization
% PrevCost = rand(ClusterNumbers,1); 
% CurrCost = zeros(ClusterNumbers,1);
% PrevPath = cell(ClusterNumbers,1);
% CurrPath = cell(ClusterNumbers,1);

SeqsNum = length(Data); 
Path = cell(SeqsNum,1);
ClusterData = cell(ClusterNumbers,1);


for I = 1:SeqsNum
    Timelength = length(Data(I).Time);
    TimeStamp = (Data(I).Stop)*IterationNum/TotalIterationNum;

    for K = 1:Timelength
        if (Data(I).Time(K) < TimeStamp)
            Data(I).Group(K) = 1;               
        else
            Data(I).Group(K) = 2;  
        end
    end
end
%Gnerating ClusterData.
 NewSeqs = BreakSeqs(Data);%break down into two seqs.
 for q = 1:ClusterNumbers
     ClusterData{q} = NewSeqs{q}; 
 end

 %10 testLocation (Exhausting Searching...)
for j = 1:ClusterNumbers
        model1(j) = Initialization_Basis(ClusterData{j});
        model(j) = Learning_MLE_S_nonstationary(ClusterData{j},model1(j),para);
end

Outputmodel = model;
ClusteredData = ClusterData;

end



