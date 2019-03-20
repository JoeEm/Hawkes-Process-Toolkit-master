%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AssginPtsToCluster() : The function receives the orignial Data , 
% then Assgins the points to corresponding cluster via Algorithnm 1,
% and finally return the ClusteredData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ClusteredData = AssginPtsToCluster(Data, ClusterNumbers, para, FirstTimeFlag)
%Initialization
PrevCost = cell(ClusterNumbers,1);
CurrCost = cell(ClusterNumbers,1);
PrevPath = cell(ClusterNumbers,1);
CurrPath = cell(ClusterNumbers,1);

%if this is the first time of the Iteration, preprocess the raw data. 
if (FirstTimeFlag == 1) 
    %e.g. if ClusterNumbers = 2, the raw data is equally divided.
    Datalength = length(Data.Time);
    if (ClusterNumbers ==2)
        for i = 1:floor((Datalength/2)) %Let first half of the Data as Cluster one 
            Data.Feature(i) = 1;
        end
        for j = (floor((Datalength/2))+1):Datalength %Let Another half of the Data as Cluster two 
            Data.Feature(i) = 2;
        end
    end                    
end

length
for i = 
    
    %Lolike_Basis()<0
     Loglike_Basis( Seqs1(1:100), model_MLE_SGL(i), para );
end
