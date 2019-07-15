%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loss() : Computing the loss of objective function 
% Rearranged Data in each cluster may not be in right order, but what we
% care about is the Time-divided boundary
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SumLoss = Loss(PointsLoss, Beta, FinalPath)
%Initialization
SumLoss = 0;
for i = 1:(length(FinalPath) - 1)
    if (FinalPath(i) == FinalPath(i + 1))
        tmp = 0;
    else
        tmp = Beta;
    end   
    SumLoss = SumLoss + tmp;
end

SumLoss = SumLoss + PointsLoss; 
end
