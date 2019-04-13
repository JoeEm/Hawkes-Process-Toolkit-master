%if the procedure doesn't stop, combine the CluterData again for convenience. 
function NewSeqs = CombineTheOriginalData(Seqs, ClusterNumbers)

ClusteredData = Seqs;
Totallength = 0;
    for g = 1:ClusterNumbers
       tmp = length(ClusteredData(g).Time);
       Totallength = tmp + Totallength;
    end
NewData = struct( 'Time',    [], ...
                                  'Mark',    [], ...
                                  'Start',   [], ...
                                  'Stop',    [], ...
                                  'Feature', [],...
                                  'Number',  []); %add Number label
    counter = 1;
    for g = 1:ClusterNumbers
       for i  = 1:length(ClusteredData(g).Time)             
           NewData.Time(counter) = ClusteredData(g).Time(i);
           NewData.Mark(counter) = ClusteredData(g).Mark(i);
          
           
           counter = counter + 1; 
       end
    end
    
NewData.Start = min(NewData.Time);
NewData.Stop =  max(NewData.Time);

NewSeqs = NewData;
    
end