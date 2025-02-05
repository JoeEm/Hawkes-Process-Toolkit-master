function Seqs = ExtractingIPTVData(Seqs, i)


%Data Processing
Seqslength = length(Seqs);
j = 1;


switch i     
        case 1  %extract 1 month data
            for k = 1:Seqslength
               while (Seqs(k).Time(j) < 30*24)

                  j = j + 1;
                  
                  if (j > length(Seqs(k).Time))
                     j = j - 1;
                     break;
                  end    
               end

               Seqs(k).Time = Seqs(k).Time(1:j);
               Seqs(k).Mark = Seqs(k).Mark(1:j);
            end


        case 2  %extract 4 months data
            for k = 1:Seqslength
               while (Seqs(k).Time(j) < 4*30*24)

                  j = j + 1;
                  
                  
                  if (j > length(Seqs(k).Time))
                     j = j - 1;
                     break;
                  end  
               end

               Seqs(k).Time = Seqs(k).Time(1:j);
               Seqs(k).Mark = Seqs(k).Mark(1:j);
            end
        
        case 3  %extract 7 months data
            for k = 1:Seqslength
               while (Seqs(k).Time(j) < 7*30*24)

                  j = j + 1;
                  
                  
                  if (j > length(Seqs(k).Time))
                     j = j - 1;
                     break;
                  end  
               end
 
               Seqs(k).Time = Seqs(k).Time(1:j);
               Seqs(k).Mark = Seqs(k).Mark(1:j);
            end
            
        case 4   %extract 10 months data    
           for k = 1:Seqslength
               while (Seqs(k).Time(j) < 10*30*24)

                  j = j + 1;
                  
                  
                  if (j > length(Seqs(k).Time))
                     j = j - 1;
                     break;
                  end  
               end

               Seqs(k).Time = Seqs(k).Time(1:j);
               Seqs(k).Mark = Seqs(k).Mark(1:j);
           end
           
        otherwise
end

end