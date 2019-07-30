%BatchSizingData by 'event'
function [NewSeqs,DSeqs] = BatchSizingData(Seqs,batchsize)

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
    while(k + batchsize -1<= Seqs(j).Stop)
    %while(k + batchsize - 1 <= length(Seqs(j).Mark))
        TimeIndex = find( (Seqs(j).Time < k + batchsize -1) & (Seqs(j).Time> k));
        TheStartIndex = TimeIndex(1);
        TheEndIndex = TimeIndex(length(TimeIndex));
        if (TheEndIndex > length(Seqs(j).Time))
            TheEndIndex = length(Seqs(j).Time);
        end
        BackupTemp(k).Time = Seqs(j).Time(TheStartIndex:TheEndIndex);
        BackupTemp(k).Mark = Seqs(j).Mark(TheStartIndex:TheEndIndex);
        BackupTemp(k).Start = min(BackupTemp(k).Time);
        BackupTemp(k).Stop = max(BackupTemp(k).Time);
        BackupTemp(k).Number = k;        
      
        k = k + 1;
    end
    NewSeqs{j} = BackupTemp;
    
%                 BackupTemp1 = struct('Time', [], ...
%       'Mark', [], ...
%       'Start', [], ...
%       'Stop', [], ...
%       'Feature', [],...
%       'Group',[]); %add Number label
    %generating New DataSeqs
%     for q = 1:minEventLength
% 
%   
%         BackupTemp1.Time(q) = Seqs(j).Time(q);
%         BackupTemp1.Mark(q) = Seqs(j).Mark(q);
%         BackupTemp1.Group(q) = Seqs(j).Group(q);  
%         
%     end
%     BackupTemp1.Start = min(BackupTemp1.Time);
%     BackupTemp1.Stop = max(BackupTemp1.Time);
%     
%     DSeqs(j) = BackupTemp1;
end
    DSeqs = Seqs;
end