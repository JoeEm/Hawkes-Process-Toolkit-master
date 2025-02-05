clear;
% 
options.N = 1000; % the number of sequences
options.Nmax = 500; % the maximum number of events per sequence
options.Tmax = 200; % the maximum size of time window
options.tstep = 0.1;
options.dt = 0.1; % the length of each time step
options.M = 250; % the number of steps in the time interval for computing sup-intensity
options.GenerationNum = 100; % the number of generations for branch processing
D = 4; % the dimension of Hawkes processes

ClusterNumbers = 3;
A1 = 0;
flag = 0;
% for j = 1:ClusterNumbers
%    %amplifier = 1 + 9*(j-1);
%     amplifier = 1;
%     MaskValue = 0.7;
% %     if (j == 3)
% %        A1 = paraDataGenrating{1}.A;
% %        flag = 1;
% %     end
% %     if (j == 4)
% %        A1 = paraDataGenrating{2}.A;
% %        flag = 1;
% %     end
%     [Seqs1{j},paraDataGenrating{j}] = GenerateSingleGroupData(options,D,j,amplifier,MaskValue,A1, flag);
%     flag = 0;
% end
% 
% 
load 730TestData.mat
InputSeqs = Seqs;


FirstTimeFlag = 1;


for k = 1:size(InputSeqs,2)
    CluNum = 1;
    Quitflag = 0;
    q=1;
    while (q <= (length(InputSeqs(k).Time)-1) )
        
      NewSeqs1{CluNum}= struct(   'Time', [], ...
                                  'Mark',    [], ...
                                  'Start',   [], ...
                                  'Stop',    [], ...
                                  'Feature', [], ...
                                  'SeqsCluNum', [], ...
                                  'Group', []);
  
       setStartpiontFlag = 1;
       setlengthSeqsFlag = 1;
       StartofaSeqsFlag = 1;
       if ( (InputSeqs(k).Group(q) ~= InputSeqs(k).Group(q+1)) && StartofaSeqsFlag == 1)
           startpoint = q;
           lengthofSeqs = 1;
           q = q+1;
           StartofaSeqsFlag = 0;
       else
           while (InputSeqs(k).Group(q) == InputSeqs(k).Group(q+1) || lengthofSeqs ==0)                  
               if (setlengthSeqsFlag == 1) 
                    lengthofSeqs = 1;
                    setlengthSeqsFlag = 0;
               end
               if (lengthofSeqs >= 1 )
                    LengthExistFlag = 0;
               end
               lengthofSeqs=lengthofSeqs+1;      
               q = q+1;
               if (q > (length(InputSeqs(k).Time)-1) )
                   Quitflag = 1;
                   break; 
               end
               if (setStartpiontFlag == 1) 
                   startpoint = q;
                   setStartpiontFlag = 0;
               end
           end
        end
%        if (LengthExistFlag == 1)
%            lengthofSeqs = 1;
%        end
       if (StartofaSeqsFlag ~= 0)
            lengthofSeqs = lengthofSeqs - 1;
       end


        counter = 1;
        for i = startpoint:( startpoint + lengthofSeqs - 1)
            NewSeqs1{CluNum}.Time(counter) = InputSeqs(k).Time(i);
            NewSeqs1{CluNum}.Mark(counter) = InputSeqs(k).Mark(i);
            NewSeqs1{CluNum}.Group(counter) = InputSeqs(k).Group(i);
            counter = counter + 1;
        end
            NewSeqs1{CluNum}.Start = min(NewSeqs1{CluNum}.Time);
            NewSeqs1{CluNum}.Stop =  max(NewSeqs1{CluNum}.Time);
            NewSeqs1{CluNum}.SeqsCluNum =  InputSeqs(k).Group(startpoint);

        if (FirstTimeFlag == 1)
               
            FHSeqs{CluNum}= struct(     'Time',    [], ...
                          'Mark',    [], ...
                          'Start',   [], ...
                          'Stop',    [], ...
                          'Feature', [], ...
                          'SeqsCluNum', [], ...
                          'Group', []);    


            for a = 1:length(NewSeqs1{CluNum}.Time)
                FHSeqs{CluNum}(k).Time(a) = NewSeqs1{CluNum}.Time(a);
                FHSeqs{CluNum}(k).Mark(a) = NewSeqs1{CluNum}.Mark(a);
                FHSeqs{CluNum}(k).Group(a) = NewSeqs1{CluNum}.Group(a);
            end
                FHSeqs{CluNum}(k).Start = min(FHSeqs{CluNum}(k).Time);
                FHSeqs{CluNum}(k).Stop =  max(FHSeqs{CluNum}(k).Time);
                FHSeqs{CluNum}(k).SeqsCluNum = NewSeqs1{CluNum}.Group(1);
        else
            if (CluNum > length( FHSeqs))
                for p  = 1:length(NewSeqs1{CluNum})
                    for o = 1:length(NewSeqs1{CluNum}(p).Time)
                        NewSeqs1{CluNum}(p).Group(o) =  FHSeqs{CluNum -1}(1).Group(1);
                    end
                end
                CluNum = CluNum - 1;               
                FHSeqs{CluNum} = [FHSeqs{CluNum},NewSeqs1{CluNum}];
            else   
                FHSeqs{CluNum} = [FHSeqs{CluNum},NewSeqs1{CluNum}];
            end
        end
        CluNum = CluNum+1;
        q = q+1;
        if (Quitflag == 1)
          break; 
        end
    end
    FirstTimeFlag = 0; 
end
    New_Seqs = FHSeqs;






