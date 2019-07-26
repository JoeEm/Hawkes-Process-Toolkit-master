function Seqs = KDDLinkingTwoSeqsToOne(Seqs1,Seqs2)
    %defalult Seqs1 & Seqs2 have the same setting

    for i = 1:length(Seqs1)
       distance = Seqs2(i).Start- Seqs1(i).Stop;
       Seqs2(i).Time = Seqs2(i).Time - distance;
       
       Seqs1(i).Time = [(Seqs1(i).Time)';  (Seqs2(i).Time)']';
       Seqs1(i).Mark = [(Seqs1(i).Mark)';  (Seqs2(i).Mark)']';
       Seqs1(i).Group = [(Seqs1(i).Group)';  (Seqs2(i).Group)']';
       
       Seqs1(i).Start = min(Seqs1(i).Time);
       Seqs1(i).Stop = max(Seqs1(i).Time);
       
    end
    
    Seqs = Seqs1;
end