%from Row52221 to Row1
%time span: 2016/6/16 : 2018/1/14
%1 week as a seq:  (7+12+1)*4 = 80 seqs.

% t = (year-1)*365*24*60 + (Month(month)+day-1)*24*60
% t = (year-1)*365 + (Month(month)+day-1);
% t = (t-time_offset)/time_scale;

load('MiniBlogData.mat');
Seqslength = length(Seqs.Time);




% for k = 1:Seqslength
%     Seqs.Time(k) = ;
% end


MiniBlogSeqs = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Feature', []);
%reverse the Seqs
for k = 1:Seqslength
    temp1 = Seqs.Time(k);
    MiniBlogSeqs.Time(k) = Seqs.Time(Seqslength - k + 1);
    MiniBlogSeqs.Time(Seqslength - k + 1) = temp1;
    
    temp2 = Seqs.Mark{k}(1);
    MiniBlogSeqs.Mark(k) = Seqs.Mark{Seqslength - k + 1}(1);
    MiniBlogSeqs.Mark(Seqslength - k + 1) = temp2;
end
MiniBlogSeqs.Time = MiniBlogSeqs.Time'; 
MiniBlogSeqs.Mark = MiniBlogSeqs.Mark'; 

year = 2016;
month = 6;
day = 16;
time_offset = MiniBlogSeqs.Time(1) - 1;

for k = 1:Seqslength
   MiniBlogSeqs.Time(k) = MiniBlogSeqs.Time(k) - time_offset;
   MiniBlogSeqs.Mark(k) = MiniBlogSeqs.Mark(k) + 1; 
end



%MiniBlogSeqsTwoYears