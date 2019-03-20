%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TotalLoss() : Computing the loss of objective function 
% Rearranged Data in each cluster may not be in right order, but what we
% care about is the Time-divided boundary
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Loss = TotalLoss( Alpha, Beta, HawkesLoglikelihood, ClusterNumbers, ClusteredData, Timeinterval)
%Initialization
Loss = 0; 
Norm2 = 0;
Tconsistency = 0;
for i = 1:ClusterNumbers
    Norm2 = Norm2 + Alpha*var(ClusteredData{i}.Time(:));
    
    Clusterlength = length(ClusteredData{i}.Time);
    tmp = sort(ClusteredData{i}.Time);
    for j = 1:Clusterlength
       a = tmp(j+1) - tmp(j);
       if (abs(a) > Timeinterval) 
           Tconsistency = Tconsistency + Beta;
       end
    end
    
    Loss = Loss + Norm2 + Tconsistency;
end

Loss = Loss + HawkesLoglikelihood; 

end