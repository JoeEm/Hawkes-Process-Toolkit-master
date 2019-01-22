%Processing IPTVData or MiniBlogRealData
%load IPTVData.mat

% Seqs1 = Seqs;
% Seqs1 = Seqs(1:30);
% Seqs1 = [Seqs1, Seqs(54:91)];
% Seqs1 = [Seqs1, Seqs(137:173)];

load MiniBlogSeqsTwoAndHalfYears.mat
Seqs1 = MiniBlogSeqsTwoAndHalfYears;

%NewSeqs = IPTV_RealData_Cutting(Seqs1);
NewSeqs = MiniBlogRealWorldCutting(Seqs1);
index = 0;
for count = 1:length(NewSeqs)
    if (sum(NewSeqs(count).Time(:)) == 0)
        index = [index, count];
    end
end

OutPutSeqs = RearrangeTheMiniBlogSeqs(NewSeqs, 4);

%TestSeqs
index2 = 1;
for count = 1:length(OutPutSeqs)
    if (sum(OutPutSeqs(count).Time(:)) == 0)
        index2 = [index2, count+1];
    end 
end

