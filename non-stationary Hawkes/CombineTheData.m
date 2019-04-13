%if the procedure doesn't stop, combine the CluterData again for convenience. 
function NewSeqs = CombineTheData(Seqs, ClusterNumbers)

NewData = cell(1,1);
for counter = 1:ClusterNumbers
    for i = 1:length(Seqs{counter})
        NewData{1} = [NewData{1},Seqs{counter}(i)] ;
    end
end


NewSeqs = NewData;
    
end