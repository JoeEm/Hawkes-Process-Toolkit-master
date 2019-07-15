function [ClusteredData,Loss, location, FinalPath] = FasterAssginPtsToCluster(Data, NewData, ClusterNumbers, para, IterationNum, Beta)

SeqsNum = length(NewData); 
ClusterData = cell(ClusterNumbers,1);

if (IterationNum == 1) %Initialization
   %computing the timeLocation for originalData
   for i = 1:ClusterNumbers - 1
       location(i) = i/3;
       TimeLocation(i) =location(i) * max(Data(1).Time);
   end
   %Partitioning the original Data(Seqs)
   flag = 0;
   for I = 1:length(Data) 
    Timelength = length(Data(I).Time);
    Q = 1;
        for K = 1:Timelength
            if (flag ~= 1)
                if ( Data(I).Time(K) < TimeLocation(Q) ) %
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
   NeSeqs = BreakSeqs(Data,ClusterNumbers);%break down into two seqs.
   for q = 1:ClusterNumbers
        ClusterData{q} = NeSeqs{q}; 
   end
   %learning the corresponding model
   for i = 1:ClusterNumbers
       model1(i) = Initialization_Basis(ClusterData{i});
       model(i) = Learning_MLE_S_nonstationary(ClusterData{i},model1(i),para);
   end
else
    %compute the new Timelocation for originalData(To be verified)
    tmplocation = cell(1,1); %first element is CluNum , second is timestamp. 
    currpos = 1;
    for j = 1:(length(NewData{1})-1) 
        if (j == 1)
            tmplocation{currpos}(1) = NewData{1}(j).ClusterNum;
            tmplocation{currpos}(2) = NewData{1}(j).Time(1);
            currpos = currpos + 1;
        end    
        if (NewData{1}(j).ClusterNum == NewData{1}(j+1).ClusterNum)
            continue;
        else
            tmplocation{currpos}(1) = NewData{1}(j+1).ClusterNum;
            tmplocation{currpos}(2) = NewData{1}(j).Time(1);
            currpos = currpos + 1;
        end
    end
    for I = 1:length(Data)
        loactionNum = length(tmplocation);
        counter = 1;
        innnerFlag = 0;
        for j = 2:loactionNum
            while( (Data(I).Time(counter) < tmplocation{j}(2)))
                Data(I).Group(counter) = tmplocation{j-1}(2);
                counter = counter + 1;
                if (counter == length(Data(I).Time))
                    break;
                    innnerFlag = 1;
                end
            end
            if (innnerFlag == 1)            
                break;
            else               
                continue;
            end
        end
    end
    
    %Partitioning the original Data(Seqs)
    flag = 0;
    for I = 1:length(Data) 
        Timelength = length(Data(I).Time);
        Q = 1;
            for K = 1:Timelength
                if (flag ~= 1)
                    if ( Data(I).Time(K) < TimeLocation(Q) ) 
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
%     NewSeqs = BreakSeqs(Data,ClusterNumbers);%break down into two seqs.
%     for q = 1:ClusterNumbers
%         ClusterData{q} = NewSeqs{q};
%     end
    
    %pcik out the longest Seqs for corresponding CluNums
    index = cell(ClusterNumbers,1);
    for i = 1:ClusterNumbers
        for j = 1:length(NewSeqs)
            if(i == NewSeqs{j}.SeqsCluNum)
              index{i} = [index{i}; i];  
            end
        end
    end
    for i = 1:ClusterNumbers
       maxIndex = 1;
       for j = 1:length(index{i})
           if ( length(NewSeqs{index{i}(j)}(1).Mark) > length(NewSeqs{index{i}(maxIndex)}(1).Mark))
                maxIndex = j;
           end
       end
       ClusterData{i} = NewSeqs{maxIndex};
    end
    
    %learning the corresponding model
    for i = 1:ClusterNumbers
        model1(i) = Initialization_Basis(ClusterData{i});
        model(i) = Learning_MLE_S_nonstationary(ClusterData{i},model1(i),para);
    end
end

%Assigning the poins to clusters.(KDD-TICC-Algorithnm 1)
PrevCost = zeros(1,ClusterNumbers);
CurrCost = zeros(1,ClusterNumbers);
PrevPath = cell(1,ClusterNumbers);
CurrPath = cell(1,ClusterNumbers);

for i = 1:length(NewData{1})
    for K = 1:ClusterNumbers
        MinIndex = find(PrevCost == min(PrevCost));
        MinIndex = MinIndex(1);
        Loglike = 0;
        for q = 1:1
            Loglike = Loglike + Loglike_Basis_NonStationary(NewData{q}(i), model(K), para);
        end
        if (PrevCost(MinIndex)) + Beta > PrevCost(K)
            CurrCost(K) = PrevCost(K) - Loglike;
            CurrPath{K} = [PrevPath{K}; K];
        else
            CurrCost(K) = PrevCost(K) - Loglike + Beta;
            CurrPath{K} = [PrevPath{K}; MinIndex];
        end
    PrevCost = CurrCost;   
    end
    PrevPath = CurrPath; 
end

FinalMinIndex = find(CurrCost == min(CurrCost));
FinalPath = CurrPath{FinalMinIndex};
Loss = CurrCost(FinalMinIndex);

%Rearrange the NewData.
for i = 1:1
   Counter = 1;
   for j = 1:length(NewData{1}) 
       NewData{i}(j).ClusterNum = FinalPath(Counter);
       Counter = Counter + 1;
   end
end


ClusteredData = NewData;

end