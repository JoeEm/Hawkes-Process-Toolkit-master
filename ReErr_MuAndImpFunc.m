%guass kernel is set as a default function.

function [RelErrMu, RelErrA] = ReErr_MuAndImpFunc( Tmax, model, modelEst)
    D = size(model.A,1);
    A = model.A;
    Aest = model.A;
    %calculating the relative error of Mu
    Mu = model.mu;
    Muest = modelEst.mu;
    w = model.w;
    RelErrMu = norm(Muest - Mu)/norm(Mu);   %The relative error of Mu
    
    ImpOptions.landmark = 0:20:100;
    ImpOptions.Tmax = Tmax;
    
    A = ImpactFunc( model, ImpOptions );
    Aest = ImpactFunc( modelEst, ImpOptions );
    
    RelErrA = norm(A(:) - Aest(:))/norm(A(:));
    
%     %calculating the relative error of ImpactFunction
%     landmark = model.landmark; %the central locations of kernels
%     sumbottom = 0;
%     sumtop = 0;
%     RelErrImFunc = 0;
%     for a = 1:D
%         for b = 1:D
%             for i = 1:Tmax
%                 distance = repmat(i, [1, length(landmark(:))]) - ...
%                 repmat(landmark(:)', [length(i), 1]);
%                 distance = distance';
%             
%                 %dbstop in ReErr_MuAndImpFunc at 23;
% 
%                 g = A(a,:,b)*(exp(-(distance.^2)./(2*w^2))./(sqrt(2*pi)*w));
%                 gest = Aest(a,:,b)*(exp(-(distance.^2)./(2*w^2))./(sqrt(2*pi)*w));
%                 temp1 = norm((gest - g),1); %To be determined
%                 
%                 sumtop = sumtop + temp1;
%                 sumbottom = sumbottom + g;
%             end
%             RelErrImFunc = RelErrImFunc + sumtop/sumbottom;    
%         end
%     end
%     RelErrImFunc = RelErrImFunc/(D^2);  %The relative error of ImpactFunction
end