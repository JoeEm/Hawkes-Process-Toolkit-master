%Processing the RealWorld Data, Let each Seqs represents a week long data.
function NewSeqs = IPTV_RealData_Cutting(Seqs)

i = 1;
j = 1;

MaxSeqsLength = length(Seqs); %1

%For IPTVData.mat only Seqs No.301 is not longer than 4 months
% index = [];

% Tmin = Seqs(1).Time(length(Seqs(1).Time - 1));
% for i = 2:MaxSeqsLength
%     t =  Seqs(i).Time(length(Seqs(i).Time) - 1);
%     if (t < Tmin)
%         Tmin = t;
% 
%     end
%      if ( t < 2880)
%          index = [index, i];
%      end
% end 
flag = 0;
HistoryEnd_i = 1;
start = 1;

NewSeqs = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Number', [],...
              'Feature', []);

j = 1;
i = 1;
Global_Count = 1;
Global_Seqs_NUM = 1;

flag  =  0; % flag of exiting the loop.
if (Seqs(1).Time(length(Seqs(1).Time)) - ...
            Seqs(1).Time(1) > 7*24)
    while (Global_Count <= MaxSeqsLength)
        j = 1;
        start = 1;
        HistoryEnd_i = 1;
        flag = 0;



        while(Seqs(Global_Count).Time(HistoryEnd_i) <= 30*10*24  || ...
                                HistoryEnd_i <= length(Seqs(Global_Count).Time))
            j = HistoryEnd_i;
            start = HistoryEnd_i;

            while ( (Seqs(Global_Count).Time(j) - Seqs(Global_Count).Time(start))  <= 7*24)
               j = j + 1; 
               if (length(Seqs(Global_Count).Time) < j )
                   fprintf(' Out of Loop\n');
                  flag  =  1;
                  break;
               end
            end

            if flag  == 0 
                NewSeqs(Global_Seqs_NUM).Time = Seqs(Global_Count).Time(start:j);
                NewSeqs(Global_Seqs_NUM).Mark = Seqs(Global_Count).Mark(start:j);
                NewSeqs(Global_Seqs_NUM).Feature = Seqs(Global_Count).Feature(start:j);
                NewSeqs(Global_Seqs_NUM).Number = Seqs(Global_Count).Number(start:j);
                NewSeqs(Global_Seqs_NUM).Start = Seqs(Global_Count).Time(start);
                NewSeqs(Global_Seqs_NUM).Stop = Seqs(Global_Count).Time(j);


                fprintf(' %d weeks sequences built already\n', Global_Seqs_NUM);
                Global_Seqs_NUM = Global_Seqs_NUM + 1;
            end
            if ( flag == 1 || j+1 >= length(Seqs(Global_Count).Time) )
                break;
            end
            HistoryEnd_i = j+1;
        end

        NewSeqs(Global_Seqs_NUM).Time = 0;  %As a interval .
        NewSeqs(Global_Seqs_NUM).Mark = 0;
        NewSeqs(Global_Seqs_NUM).Start = 0;
        NewSeqs(Global_Seqs_NUM).Stop = 0;
        Global_Seqs_NUM = Global_Seqs_NUM + 1;
        Global_Count = Global_Count + 1;
        %dbstop in IPTV_RealData_Cutting at 78;
        fprintf('broke down %d sequences already\n', Global_Count);
    end
else
    NewSeqs = Seqs;
end
    
end