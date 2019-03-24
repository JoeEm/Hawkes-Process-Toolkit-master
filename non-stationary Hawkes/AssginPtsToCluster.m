%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AssginPtsToCluster() : The function receives the orignial Data , 
% then Assgins the points to corresponding cluster via Algorithnm 1,
% and finally return the ClusteredData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ClusteredData = AssginPtsToCluster(Data, ClusterNumbers, para, FirstTimeFlag, Beta)
%Initialization
PrevCost = rand(ClusterNumbers,1); 
CurrCost = zeros(ClusterNumbers,1);
DataBatchLength = 50;


%if this is the first time of the Iteration, preprocess the raw data. 
if (FirstTimeFlag == 1) 
    %e.g. if ClusterNumbers = 2, the raw data is equally divided.
    Datalength = length(Data.Time);
    if (ClusterNumbers ==2)
        for i = 1:floor((Datalength/2)) %Let first half of the Data as Cluster one 
            Data.Feature(i) = 1;
        end
        for j = (floor((Datalength/2))+1):Datalength %Let Another half of the Data as Cluster two 
            Data.Feature(j) = 2;
        end
    end                    
end



start = zeros(ClusterNumbers,1);
stop  = zeros(ClusterNumbers,1);
Path = zeros(Datalength - DataBatchLength,1) + 1;
PrevSelection = 1;

%for ClusterNum = 2 situation
BoundaryIndex = 0;
for i = 1:Datalength
    if Data.Feature(i) == 1
       BoundaryIndex = BoundaryIndex + 1;
    end
end
start(1,1) = 1;
start(2,BoundaryIndex) = BoundaryIndex;
stop(1,1) = BoundaryIndex+1;
stop(2,BoundaryIndex) = Datalength;

%dbstop in AssginPtsToCluster at 46
DataBackup = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Feature', []);
          
for i = 1:(Datalength - DataBatchLength)
    for j = 1:ClusterNumbers
        MinIndex = find(PrevCost == min(PrevCost)); %find the smallest index
        %
        counter = 1; 
        for j = i:(i + DataBatchLength)
            DataBackup.Time(counter) = Data.Time(j);
            DataBackup.Mark(counter) = Data.Mark(j);
            DataBackup.Feature(counter) = Data.Feature(j);
            counter = counter + 1;
        end
        %
        model = Initialization_Basis(Data);
        model = Learning_MLE_S_Basis(Data , model, para);
        Logvalue = Loglike_Basis_NonStationary( DataBackup, model, para, start, stop);
        if (PrevCost(MinIndex) + Beta > PrevCost(j))
            CurrCost(j) = PrevCost(j) - Logvalue + Beta;
            Path(i) = j;
            PrevSelection = j;
        else
            CurrCost(j) = PrevCost(MinIndex) - Logvalue;
            Path(i) = PrevSelection;
        end   
        PrevCost = CurrCost;
    end
end

tmep = cell(ClusterNumbers,1) ;
for i = 1:ClusterNumbers
    temp{i} = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Feature', []);
end

% specification: [array,index] = sort([2 3 1])  index = [3 1 2] array =[1 2 3]
for i = 1:(Datalength - DataBatchLength)
    ClusterNum = Path(i);
    tmp{ClusterNum}.Time = [tmp{ClusterNum}.Time, Data(i).Time];
    tmp{ClusterNum}.Mark = [tmp{ClusterNum}.Mark, Data(i).Mark];
    tmp{ClusterNum}.Feature = [tmp{ClusterNum}.Feature, Data(i).Feature];
end%the series data in each cluster is out of order.

%rearrange the Cluster
BackupTemp{ClusterNum} = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Feature', []);
for i = 1:ClusterNum
    %[array,index] = sort([2 3 1])
    [NewTime, index] = sort(tmp{ClusterNum}.Time);
    for j = 1:length(index)
        BackupTemp{ClusterNum}.Time(j) = tmp{ClusterNum}.Time(index(j));
        BackupTemp{ClusterNum}.Mark(j) = tmp{ClusterNum}.Mark(index(j));
        BackupTemp{ClusterNum}.Feature(j) = tmp{ClusterNum}.Feature(index(j));
    end
end

ClusteredData = BackupTemp;  

end

