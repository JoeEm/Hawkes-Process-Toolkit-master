%add Number label to each data for later computing of loss
function NewSeqs = NumberTheData(Seqs)
BackupTemp = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Feature', [],...
              'Number',[]); %add Number label

if (~isempty(Seqs.Feature))
    for j = 1:length(Seqs.Time)
         BackupTemp.Time(j) = Seqs.Time(j);
         BackupTemp.Mark(j) = Seqs.Mark(j);
         BackupTemp.Feature(j) = Seqs.Feature(j);
         BackupTemp.Number(j) = j;
    end
else
    for j = 1:length(Seqs.Time)
         BackupTemp.Time(j) = Seqs.Time(j);
         BackupTemp.Mark(j) = Seqs.Mark(j);

         BackupTemp.Number(j) = j;
    end
end
     BackupTemp.Start = min(BackupTemp.Time);
     BackupTemp.Stop = max(BackupTemp.Time); 
   
     NewSeqs = BackupTemp;
end