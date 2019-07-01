%testing the result-A of MCMC method.
load FirstHalfTest.mat

T = size(A,3);
D=size(A,2);
totalNUm = D*D;

% for k = 1:T
%     threshold = max(max(A(:,:,k)));
%     for i = 1:D
%         for j = 1:D
%             if ( A(i,j,k) < threshold )
%                 A(i,j,k) = 0;
%             end
%         end
%     end
% end


A_New = A(:,:,2:T);

% Visualize the infectivity matrix (the adjacent matrix of Granger causality graph)
q=1;
while(1)
    figure
    subplot(221)        
    imagesc(A_New(:,:,q))
    title(sprintf('T = %d',q));
    axis square
    colorbar

    subplot(222)         
    imagesc(A_New(:,:,q+1))
    title(sprintf('T = %d',q+1));
    colorbar
    axis square

    subplot(223)        
    imagesc(A_New(:,:,q+2))
    title(sprintf('T = %d',q+2));
    colorbar
    axis square

    subplot(224)        
    imagesc(A_New(:,:,q+3))
    title(sprintf('T = %d',q+3));
    colorbar
    axis square

    
    q = q+4;
    if (q > T)
       break; 
    end
end