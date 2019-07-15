%BatchSizingData by 'event'
function NewSeqs = BatchSizingData(Seqs,batchsize)

NewSeqs = cell(1,size(Seqs,2));
%find the Seqs with the least 'event' Numbers.
% minIndex = 1;
% minEventLength = length(Seqs(minIndex).Mark);
% for i = 1:size(Seqs,2)
%     if ( length(Seqs(i).Mark) < minEventLength)
%         minIndex = i;
%         minEventLength = length(Seqs(minIndex).Mark);
%     end
% end

for j = 1:size(Seqs,2)
    k =1;
            
    BackupTemp = struct('Time', [], ...
      'Mark', [], ...
      'Start', [], ...
      'Stop', [], ...
      'ClusterNum', [],...
      'Number',[]); %add Number label
%     while(k + batchsize - 1 <= minEventLength)
    while(k + batchsize - 1 <= length(Seqs(j).Mark))
        BackupTemp(k).Time = Seqs(j).Time(k:(k + batchsize - 1));
        BackupTemp(k).Mark = Seqs(j).Mark(k:(k + batchsize - 1));
        BackupTemp(k).Start = min(BackupTemp(k).Time);
        BackupTemp(k).Stop = max(BackupTemp(k).Time);
        BackupTemp(k).Number = k;        
      
        k = k + 1;
    end
    NewSeqs{j} = BackupTemp;
end
    
end