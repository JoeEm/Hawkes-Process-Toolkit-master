function NewSeqs = CutSeqs(InputSeqs, length)

NewSeqs= struct(     'Time',    [], ...
                      'Mark',    [], ...
                      'Start',   [], ...
                      'Stop',    [], ...
                      'Feature', []);
for i = 1:length
    NewSeqs.Time(i) = InputSeqs.Time(i);
    NewSeqs.Mark(i) = InputSeqs.Mark(i);
end
    
NewSeqs.Start = min(NewSeqs.Time);
NewSeqs.Stop =  max(NewSeqs.Time);

end