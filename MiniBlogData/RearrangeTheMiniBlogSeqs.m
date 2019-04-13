function OutPutSeqs = RearrangeTheMiniBlogSeqs(NewSeqs, i)

OutPutSeqs = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Feature', []);

Seqslength = length(NewSeqs);
%find the Separator Index
index = 1;
 fprintf(' Seqslength = %d \n', Seqslength);           
        
for k = 1:Seqslength
   if ( sum(NewSeqs(k).Mark) == 0)
       if (k ~= Seqslength)
          index = [index, k+1];
       end
   end
end

SegmentNum = length(index);
save('index');
 fprintf(' SegmentNum = %d \n', SegmentNum);

if (length(index) > 1)
    %rearrange the Seqs    
    interval = 1;
    start = index(1);
    OutPutSeqs = NewSeqs(start:index(2)-2);

        for k = 2:SegmentNum
            start = index(k);
            fprintf(' start = %d \n', start);
            %dbstop in RearrangeTheSeqs at 35;
    %         if (start + interval  > length(NewSeqs))
    %            break; 
    %         end
            if k == SegmentNum
                OutPutSeqs = [OutPutSeqs, NewSeqs( start: end-1 )];  
            else   
                OutPutSeqs = [OutPutSeqs, NewSeqs( start: index(k+1)-2 )];
            end
        end

    %     %Rearrange the Data again
    %     index2 = 1;
    %     for count = 1:length(OutPutSeqs)
    %         if (sum(OutPutSeqs(count).Time(:)) == 0)
    %             index2 = [index2, count+1];
    %         end 
    %     end
    %     
    %     for k = 2:length(index2)
    %         location = index2(k);
    %         if k ~= length(index2)
    %             for g = location: index2(k+1) - 2
    %                 OutPutSeqs(g - 1) = OutPutSeqs(g);
    %                 g = g+1;
    %             end
    %         else
    %             for g = location: length(OutPutSeqs)-1
    %                 OutPutSeqs(g - 1) = OutPutSeqs(g);
    %                 g = g+1;
    %             end
    %         end
    %         
    %         
    %     end
else
    if (length(NewSeqs) > 1)
        OutPutSeqs = NewSeqs(1:length(NewSeqs)-1);
    else
        OutPutSeqs = NewSeqs;
    end
end




end

