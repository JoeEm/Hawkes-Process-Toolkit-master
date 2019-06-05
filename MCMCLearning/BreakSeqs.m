function [FHSeqs1, FHSeqs2]  = BreakSeqs(InputSeqs)

FirstTimeFlag = 1;

for k = 1:size(InputSeqs,2)


        NewSeqs1= struct(     'Time',    [], ...
                          'Mark',    [], ...
                          'Start',   [], ...
                          'Stop',    [], ...
                          'Feature', []);
    
    lengthofSeqs = sum(InputSeqs(k).Time <= 200);
    
    for i = 1:lengthofSeqs
        NewSeqs1.Time(i) = InputSeqs(k).Time(i);
        NewSeqs1.Mark(i) = InputSeqs(k).Mark(i);
    end

    NewSeqs1.Start = min(NewSeqs1.Time);
    NewSeqs1.Stop =  max(NewSeqs1.Time);


    NewSeqs2= struct(     'Time',    [], ...
                          'Mark',    [], ...
                          'Start',   [], ...
                          'Stop',    [], ...
                          'Feature', []);
                      counter = 1;
    for j = lengthofSeqs+1:length(InputSeqs(k).Time)
        NewSeqs2.Time(counter) = InputSeqs(k).Time(j);
        NewSeqs2.Mark(counter) = InputSeqs(k).Mark(j);
        counter = counter + 1 ;
    end

    NewSeqs2.Start = min(NewSeqs2.Time);
    NewSeqs2.Stop =  max(NewSeqs2.Time);
    
    if (FirstTimeFlag == 1)
       FHSeqs1= struct(     'Time',    [], ...
                      'Mark',    [], ...
                      'Start',   [], ...
                      'Stop',    [], ...
                      'Feature', []);
        FHSeqs2= struct(     'Time',    [], ...
                      'Mark',    [], ...
                      'Start',   [], ...
                      'Stop',    [], ...
                      'Feature', []);
        for a = 1:length(NewSeqs1.Time)
            FHSeqs1.Time(a) = NewSeqs1.Time(a);
            FHSeqs1.Mark(a) = NewSeqs1.Mark(a);
        end
        
            FHSeqs1.Start = min(FHSeqs1.Time);
            FHSeqs1.Stop =  max(FHSeqs1.Time);
        
        for b = 1:length(NewSeqs2.Time)
            FHSeqs2.Time(b) = NewSeqs2.Time(b);
            FHSeqs2.Mark(b) = NewSeqs2.Mark(b);
        end
            FHSeqs2.Start = min(FHSeqs2.Time);
            FHSeqs2.Stop =  max(FHSeqs2.Time);        
        FirstTimeFlag = 0;
    else
        FHSeqs1 = [FHSeqs1,NewSeqs1];
        FHSeqs2 = [FHSeqs2,NewSeqs2];
    end 


  
    

end
    
end