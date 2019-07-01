function BatchData  = CombinePointsToOneBatch(InputSeqs,BatchSize)

Seqslength = size(InputSeqs(1).Time);
InterationNum = ceil(Seqslength/BatchSize);

GlobalSeqslengthCounter = 1;
GlobalDataLocationCounter = 1;
EndFlag = 0;

for k = 1:InterationNum

    TmpS = struct(     'Time',    [], ...
                      'Mark',    [], ...
                      'Start',   [], ...
                      'Stop',    [], ...
                      'Feature', [], ...
                      'Group', []);
    


    counter = 1;
    while(1)
        if (GlobalSeqslengthCounter > Seqslength)
            EndFlag = 1; 
        end

        if (counter > BatachSize || EndFlag == 1)
           break; 
        end
        
        TmpS.Time(counter) = InputSeqs.Time(GlobalDataLocationCounter);
        TmpS.Mark(counter) = InputSeqs.Mark(GlobalDataLocationCounter);
        TmpS.Group(counter) = InputSeqs.Group(GlobalDataLocationCounter);
        counter = counter + 1;
        GlobalDataLocationCounter = GlobalDataLocationCounter + 1;
        GlobalSeqslengthCounter = GlobalSeqslengthCounter + 1;
    end

    if (EndFlag == 1)
        break; 
    end
    TmpS.Start = min(TmpS.Time);
    TmpS.Stop =  max(TmpS.Time);

    BatchData(k) =  TmpS;
    
    TmpS = [];%clear the structure.

end