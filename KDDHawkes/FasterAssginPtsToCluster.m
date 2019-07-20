function [ClusteredData,Loss, tmplocation, FinalPath,Data] = FasterAssginPtsToCluster(Data, NewData, ClusterNumbers, para, IterationNum, Beta, batchsize)

SeqsNum = length(NewData); 
ClusterData = cell(ClusterNumbers,1);

if (IterationNum == 1) %Initialization
   %computing the timeLocation for originalData
   for i = 1:ClusterNumbers - 1
       location(i) = i/3;
       TimeLocation(i) =location(i) * max(Data(1).Time);
   end
   tmplocation = TimeLocation;
   %Partitioning the original Data(Seqs)
   flag = 0;
   for I = 1:length(Data) 
    Timelength = length(Data(I).Time);
    Q = 1;
        for K = 1:Timelength
            if (flag ~= 1)
                if( K < (location(Q)*length(Data(I).Mark)))
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
   NeSeqs = BreakSeqs(Data);%break down into two seqs.
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
            tmplocation{currpos}(2) = NewData{1}(j).Time(1+batchsize/2);
            currpos = currpos + 1;
        end
        if (j + 1  == length(length(NewData{1})))
            currpos = currpos + 1;
            tmplocation{currpos}(1) = NewData{1}(j).ClusterNum;
            tmplocation{currpos}(2) = NewData{1}(j).Time(1+batchsize/2);
            continue;
        end
        if (NewData{1}(j).ClusterNum == NewData{1}(j+1).ClusterNum)
            continue;
        else
            tmplocation{currpos}(1) = NewData{1}(j+1).ClusterNum;
            tmplocation{currpos}(2) = NewData{1}(j).Time(1+batchsize/2);
            currpos = currpos + 1;
        end
    end
    for I = 1:length(Data)
        loactionNum = length(tmplocation);
        counter = 1;
        innnerFlag = 0;
        for j = 2:loactionNum
            while( counter < ((tmplocation{j}(2))/((Data(1).Stop-Data(1).Start))*length(Data(1).Mark)))
                Data(I).Group(counter) = tmplocation{j-1}(1);
                counter = counter + 1;
                if (counter == length(Data(I).Time))
                    innnerFlag = 1;
                    break;
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
%     flag = 0;
%     for I = 1:length(Data) 
%         Timelength = length(Data(I).Time);
%         Q = 1;
%             for K = 1:Timelength
%                 if (flag ~= 1)
%                     if (  K < ((tmplocation{Q}(2))/((Data(I).Stop-Data(I).Start))*length(Data(I).Mark)) ) 
%                         Data(I).Group(K) = tmplocation{Q}(1);                               
%                     else
%                         if (Q+1 == length(tmplocation))
%                             flag = 1;
%                             Q = Q+1;
%                         else                         
%                             Q = Q+1;
%                         end
%                     end
%                 else
%                     Data(I).Group(K) = tmplocation{Q}(1);    
%                 end
%             end
%             flag = 0;
%     end
    
    %Gnerating ClusterData.
    NewSeqs = BreakSeqs(Data);%break down into two seqs.
%     for q = 1:ClusterNumbers
%         ClusterData{q} = NewSeqs{q};
%     end
    
    %pcik out the longest Seqs for corresponding CluNums
    index = cell(ClusterNumbers,1);
    for i = 1:ClusterNumbers
        for j = 1:length(NewSeqs)
            if(i == NewSeqs{j}(1).SeqsCluNum)
              index{i} = [index{i}; j];  
            end
        end
    end
    for i = 1:ClusterNumbers
       maxIndex = index{i}(1);
       for j = 1:length(index{i})
           if ( length(NewSeqs{index{i}(j)}(1).Mark) > length(NewSeqs{maxIndex}(1).Mark))
                maxIndex = index{i}(j);
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

%Assigning the poins to clusters.(KDD-TICC-Algorithnm 1) greedy algorithnm
Path = [];
CurrCost = 0;
for i = 1:length(NewData{1})
    minindex = 1;
    Loglike(1,ClusterNumbers) = 0;
    for K = 1:ClusterNumbers
        Loglike(1,K) = -Loglike_Basis_NonStationary(NewData{1}(i), model(K), para);
    end
    index = find(Loglike == min(Loglike));
    CurrCost = CurrCost + Loglike(index(1));
    Path = [Path;index];
end

FinalPath = Path;
Loss = CurrCost;

% 
%Assigning the poins to clusters.(KDD-TICC-Algorithnm 1) dynamic programming algorithnm
% PrevCost = zeros(1,ClusterNumbers);
% CurrCost = zeros(1,ClusterNumbers);
% PrevPath = cell(1,ClusterNumbers);
% CurrPath = cell(1,ClusterNumbers);
%tmp = zeros(ClusterNumbers,length(NewData{1}));

% for i = 1:length(NewData{1})
%     for K = 1:ClusterNumbers
%         MinIndex = find(PrevCost == min(PrevCost));
%         MinIndex = MinIndex(1);
%         Loglike = 0;
%         for q = 1:1
%             Loglike = Loglike + Loglike_Basis_NonStationary(NewData{q}(i), model(K), para);
%             tmp(K,i) = Loglike;
%         end
%         if (PrevCost(MinIndex)) + Beta > PrevCost(K)
%             CurrCost(K) = PrevCost(K) - Loglike;
%             CurrPath{K} = [PrevPath{K}; K];
%         else
%             CurrCost(K) = PrevCost(K) - Loglike + Beta;
%             CurrPath{K} = [PrevPath{K}; MinIndex];
%         end
%     PrevCost = CurrCost;   
%     end
%     PrevPath = CurrPath; 
% end
% 

% FinalMinIndex = find(CurrCost == min(CurrCost));
% FinalPath = CurrPath{FinalMinIndex};
%Loss = CurrCost(FinalMinIndex);

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