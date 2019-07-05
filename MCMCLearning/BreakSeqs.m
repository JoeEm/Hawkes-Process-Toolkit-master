function New_Seqs  = BreakSeqs(InputSeqs,ClusterNumbers)
FirstTimeFlag = 1;

for CluNum = 1:ClusterNumbers    
    FHSeqs{CluNum}= struct(     'Time',    [], ...
                  'Mark',    [], ...
                  'Start',   [], ...
                  'Stop',    [], ...
                  'Feature', [], ...
                  'Group', []);    
end

for k = 1:size(InputSeqs,2)   
    for CluNum = 1:ClusterNumbers
        NewSeqs1{CluNum}= struct(     'Time',    [], ...
              'Mark',    [], ...
              'Start',   [], ...
              'Stop',    [], ...
              'Feature', [], ...
              'Group', []);

        if (CluNum == 1)
            lengthofSeqs = sum(InputSeqs(k).Group == CluNum);
            startpoint = 1;
        else
            startpoint = startpoint + lengthofSeqs;
            lengthofSeqs = sum(InputSeqs(k).Group == CluNum);%problem sentence
        end
        counter = 1;
        for i = startpoint:( startpoint + lengthofSeqs - 1)
            NewSeqs1{CluNum}.Time(counter) = InputSeqs(k).Time(i);
            NewSeqs1{CluNum}.Mark(counter) = InputSeqs(k).Mark(i);
            NewSeqs1{CluNum}.Group(counter) = InputSeqs(k).Group(i);
            counter = counter + 1;
        end
            NewSeqs1{CluNum}.Start = min(NewSeqs1{CluNum}.Time);
            NewSeqs1{CluNum}.Stop =  max(NewSeqs1{CluNum}.Time);
        
        if (FirstTimeFlag == 1)
            for a = 1:length(NewSeqs1{CluNum}.Time)
                FHSeqs{CluNum}(k).Time(a) = NewSeqs1{CluNum}.Time(a);
                FHSeqs{CluNum}(k).Mark(a) = NewSeqs1{CluNum}.Mark(a);
                FHSeqs{CluNum}(k).Group(a) = NewSeqs1{CluNum}.Group(a);
            end
                FHSeqs{CluNum}(k).Start = min(FHSeqs{CluNum}(k).Time);
                FHSeqs{CluNum}(k).Stop =  max(FHSeqs{CluNum}(k).Time);       
        else       
            FHSeqs{CluNum} = [FHSeqs{CluNum},NewSeqs1{CluNum}];
        end            
    end
    FirstTimeFlag = 0;
end
    New_Seqs = FHSeqs;
end