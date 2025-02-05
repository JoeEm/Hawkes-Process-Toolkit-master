%add Number label to each Seqs for later computing of loss
%initialize the 'Feture' (clusterNum that Seq belongs to)
function NewSeqs = NumberTheSeqs(Seqs)
BackupTemp = struct('Time', [], ...
          'Mark', [], ...
          'Start', [], ...
          'Stop', [], ...
          'Feature', [],...
          'Number',[]); %add Number label
TempSeqs = struct('Time', [], ...
          'Mark', [], ...
          'Start', [], ...
          'Stop', [], ...
          'Feature', [],...
          'Number',[]); %add Number label

Datalength = length(Seqs);
%add Number label to each Seqs(Initialization)
for i = 1:floor(Datalength/4) % 1-1000  -> 1
    Seqs(i).Feature = 1;
end
for j = (floor((Datalength/4))+1):floor(Datalength/2) % 1001-2000 ->2
    Seqs(j).Feature = 2;
end 
for j = (floor((Datalength/2))+1):floor(Datalength*3/4) % 2001-3000 ->2
    Seqs(j).Feature = 2;
end
for j = (floor(Datalength*3/4)+1):Datalength % 3001-4000 ->1
    Seqs(j).Feature = 2;
end 

% for i = 1:floor(Datalength/2) % 1-2000  -> 1
%     Seqs(i).Feature = 1;
% end
% 
% for j = (floor(Datalength/2)+1):Datalength % 2001-4000 ->2
%     Seqs(j).Feature = 2;
% end 
% 

for i = 1:length(Seqs)
    BackupTemp.Time = Seqs(i).Time;
    BackupTemp.Mark = Seqs(i).Mark;
    BackupTemp.Start = Seqs(i).Start;
    BackupTemp.Stop = Seqs(i).Stop;
    BackupTemp.Feature = Seqs(i).Feature;
    BackupTemp.Number = i;

    if (i == 1)
        TempSeqs = BackupTemp;
    else
        TempSeqs = [TempSeqs, BackupTemp];
    end
    BackupTemp.Time = [];
    BackupTemp.Mark = [];
    BackupTemp.Start = [];
    BackupTemp.Stop = [];
    BackupTemp.Feature = [];
    BackupTemp.Number = [];
end

                   
    NewSeqs = TempSeqs;
end
