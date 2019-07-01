function New_Seqs  = BreakSeqs(InputSeqs,options)

FirstTimeFlag = 1;

for k = 1:size(InputSeqs,2)


        NewSeqs1= struct(     'Time',    [], ...
                          'Mark',    [], ...
                          'Start',   [], ...
                          'Stop',    [], ...
                          'Feature', [], ...
                          'Group', []);
    
    lengthofSeqs = sum(InputSeqs(k).Group == 1);
    
    for i = 1:lengthofSeqs
        NewSeqs1.Time(i) = InputSeqs(k).Time(i);
        NewSeqs1.Mark(i) = InputSeqs(k).Mark(i);
        NewSeqs1.Group(i) = InputSeqs(k).Group(i);
    end

    NewSeqs1.Start = min(NewSeqs1.Time);
    NewSeqs1.Stop =  max(NewSeqs1.Time);

    NewSeqs2= struct(     'Time',    [], ...
                          'Mark',    [], ...
                          'Start',   [], ...
                          'Stop',    [], ...
                          'Feature', [], ...
                          'Group', []);
                      counter = 1;
    for j = lengthofSeqs+1:length(InputSeqs(k).Time)
        NewSeqs2.Time(counter) = InputSeqs(k).Time(j);
        NewSeqs2.Mark(counter) = InputSeqs(k).Mark(j);
        NewSeqs2.Group(counter) = InputSeqs(k).Group(j);
        counter = counter + 1 ;
    end

    NewSeqs2.Start = min(NewSeqs2.Time);
    NewSeqs2.Stop =  max(NewSeqs2.Time);
    
    if (FirstTimeFlag == 1)
       FHSeqs1= struct(     'Time',    [], ...
                      'Mark',    [], ...
                      'Start',   [], ...
                      'Stop',    [], ...
                      'Feature', [], ...
                      'Group', []);
        FHSeqs2= struct(     'Time',    [], ...
                      'Mark',    [], ...
                      'Start',   [], ...
                      'Stop',    [], ...
                      'Feature', [], ...
                      'Group', []);
        for a = 1:length(NewSeqs1.Time)
            FHSeqs1.Time(a) = NewSeqs1.Time(a);
            FHSeqs1.Mark(a) = NewSeqs1.Mark(a);
            FHSeqs1.Group(a) = NewSeqs1.Group(a);
        end
            FHSeqs1.Start = min(FHSeqs1.Time);
            FHSeqs1.Stop =  max(FHSeqs1.Time);
        
        for b = 1:length(NewSeqs2.Time)
            FHSeqs2.Time(b) = NewSeqs2.Time(b);
            FHSeqs2.Mark(b) = NewSeqs2.Mark(b);
            FHSeqs2.Group(b) = NewSeqs2.Group(b);
        end
            FHSeqs2.Start = min(FHSeqs2.Time);
            FHSeqs2.Stop =  max(FHSeqs2.Time);        
        FirstTimeFlag = 0;
    else
        FHSeqs1 = [FHSeqs1,NewSeqs1];
        FHSeqs2 = [FHSeqs2,NewSeqs2];
    end 

end
    New_Seqs{1} =  FHSeqs1;
    New_Seqs{2} =  FHSeqs2; 
    
end