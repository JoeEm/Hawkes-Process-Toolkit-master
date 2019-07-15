%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TotalLoss() : Computing the loss of objective function 
% Rearranged Data in each cluster may not be in right order, but what we
% care about is the Time-divided boundary
% ClusteredData{ClusterNumbers}  cell structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Loss,Loglike] = TotalLoss( ClusterNumbers, ClusteredData, model, alg)
%Initialization
% Loss = 0; 
% Norm2 = 0;
% Tconsistency = 0;
% for i = 1:ClusterNumbers
%     for q = 1:length(ClusteredData{i})
%        b(q) = ClusteredData{i}(q).Number;
%     end
%     %sorted the
%     [array,index] = sort(b); 
%     %Norm2 = Norm2 + Alpha*var(ClusteredData{i}.Time(:));
%     
%     Clusterlength = length(b);
%     %[array,index] = sort([2 3 1])  index = [3 1 2] array =[1 2 3]
%     %ClusterData should've been already sorted from the dimension of 'Time'.
%     for j = 1:Clusterlength - 1
%        a = array(j+1) - array(j);
%        if (abs(a) > 1) %neighbouring data point is supposed to be assgined to the same cluster.
%            Tconsistency = Tconsistency + Beta;
%        end
%     end
%     
%     Loss = Loss + Norm2 + Tconsistency;
% end
% 
% Loss = Loss + HawkesLoglikelihood; 
Loss = 0;
for i = 1:ClusterNumbers
    Loglike(i) = Loglike_Basis_NonStationary(ClusteredData{i}, model(i), alg);
    Loss = Loss + Loglike(i);
end

end
