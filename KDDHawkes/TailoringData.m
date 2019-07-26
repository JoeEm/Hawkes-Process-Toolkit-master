function Seqs1 = TailoringData(Seqs)

%find the Seqs with the least 'event' Numbers.
minIndex = 1;
minEventLength = length(Seqs(minIndex).Mark);
for i = 1:size(Seqs,2)
    if ( length(Seqs(i).Mark) < minEventLength)
        minIndex = i;
        minEventLength = length(Seqs(minIndex).Mark);
    end
end

for j = 1:size(Seqs,2)
            
    BackupTemp(j) = struct('Time', [], ...
      'Mark', [], ...
      'Start', [], ...
      'Stop', [], ...
      'Group', []); %add Number label

    %while(k + batchsize - 1 <= length(Seqs(j).Mark))
    BackupTemp(j).Time = Seqs(j).Time(1:minEventLength);
    BackupTemp(j).Mark = Seqs(j).Mark(1:minEventLength);
    BackupTemp(j).Group = Seqs(j).Group(1:minEventLength);
    BackupTemp(j).Start = min(BackupTemp(j).Time);
    BackupTemp(j).Stop = max(BackupTemp(j).Time);
    
end
    Seqs1 = BackupTemp;
end