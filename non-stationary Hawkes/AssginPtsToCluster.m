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

for Q = 1:ClusterNumbers-1 %e.g a seqs is separated into n blcoks through n-1   Tmax = 200
   location(Q) = Q/ClusterNumbers; 
   tmp = ((Data(1).Stop)/(Q+1));
   TimeStamp(Q) = (Data(1).Stop)*location(Q) + tmp*(IterationNum/TotalIterationNum);
end

flag = 0;
for I = 1:SeqsNum %Partition strategy
    Timelength = length(Data(I).Time);
    Q = 1;
    for K = 1:Timelength
        if (flag ~= 1)
            if ( Data(I).Time(K) < TimeStamp(Q) ) %
                Data(I).Group(K) = Q;                               
            else
                if (Q ~= ClusterNumbers - 1)
                    Q = Q+1;
                else
                    flag = 1;
                    Q = Q+1;
                end
            end
        else
            Data(I).Group(K) = Q;    
        end
    end
    flag = 0;
end
%Gnerating ClusterData.
 NewSeqs = BreakSeqs(Data,ClusterNumbers);%break down into two seqs.
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



