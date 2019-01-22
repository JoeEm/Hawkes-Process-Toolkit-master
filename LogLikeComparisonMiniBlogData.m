clear %clear the workplace Data

options.Tmax  = 7;

load('MiniBlogSeqsTwoAndHalfYears.mat');

NewSeqs = MiniBlogRealWorldCutting(MiniBlogSeqsTwoAndHalfYears);
OutPutSeqs = RearrangeTheMiniBlogSeqs(NewSeqs, 4);

%NewSeqs = NewSeqs(1:length(NewSeqs) - 1);
Seqs1 = OutPutSeqs;


disp('Maximum likelihood estimation and basis representation') 
alg1.LowRank = 0; % without low-rank regularizer
alg1.Sparse = 1; % with sparse regul  arizer
alg1.alphaS = 1;
alg1.GroupSparse = 1; % with group-sparse regularizer
alg1.alphaGS = 100;
alg1.outer = 8; %8
alg1.rho = 0.1; % the initial parameter for ADMM
alg1.inner = 5; %5
alg1.thres = 1e-5;
alg1.Tmax = [];
alg1.storeErr = 0;
alg1.storeLL = 0;

alg1.As = 10; %for Local Independence R.
alg1.Ap = 1000;%for Pairwise similarity R.



model1 = Initialization_Basis(Seqs1); 
%model1.kernel = 'gauss';
%model1.w = pi*800/(7);  %According to the Real-Data Experiment
%model1.landmark = (model1.w)*(0:ceil(options.Tmax/model1.w));
%D = 7 ;
%model1.A = rand(D, length(model1.landmark), D);
%model1.mu = rand(D, 1)./D;


% learning the model by MLE
alg1.Sparse=0;
alg1.GroupSparse=0; 
model = Learning_MLE_Basis( Seqs1, model1, alg1 ); 
A1 = ImpactFunc( model, options );

% learning the model by MLE-S
alg1.Sparse=1;
alg1.GroupSparse=0; 
model_S = Learning_MLE_S_Basis( Seqs1, model1, alg1 ); 
A1_S = ImpactFunc( model_S, options );

% learning the model by MLE-SGL
alg1.Sparse=1;
alg1.GroupSparse=1; 
model_MLE_SGL = Learning_MLE_S_Basis( Seqs1, model1, alg1 );
A1_SGL = ImpactFunc(model_MLE_SGL, options);


% Visualize the infectivity matrix (the adjacent matrix of Granger causality graph)
figure
subplot(131)        
imagesc(A1)
title('Estimated infectivity-MLE') 
axis square
colorbar

subplot(132)        
imagesc(A1_S)
title('Estimated infectivity-MLE-S')
colorbar
axis square

subplot(133)        
imagesc(A1_SGL)
title('Estimated infectivity-MLE-SGL')
colorbar
axis square




% nn = 1; %5iteration
% for i = 1:nn
%     base = 500;
%     TrainNum = base*i;
% 
%     % learning the model by MLE
%     alg1.Sparse=0;
%     alg1.GroupSparse=0; 
%     model(i) = Learning_MLE_Basis( Seqs1(1:TrainNum), model1, alg1 ); 
%     A1{i} = ImpactFunc( model(i), options );
%     LogLike_MLE(i) = Loglike_Basis( Seqs1(2001:3000), model(i), alg1 );
% 
% 
%     % learning the model by MLE-S
%     alg1.Sparse=1;
%     alg1.GroupSparse=0; 
%     model_S(i) = Learning_MLE_S_Basis( Seqs1(1:TrainNum), model1, alg1 ); 
%     A1_S{i} = ImpactFunc( model_S(i), options );
%     LogLike_MLE_S(i) = Loglike_Basis( Seqs1(2001:3000), model_S(i), alg1 );
% 
% 
%     % learning the model by MLE-SGL
%     alg1.Sparse=1;
%     alg1.GroupSparse=1; 
%     model_SGL(i) = Learning_MLE_S_Basis( Seqs1(1:TrainNum), model1, alg1 ); 
%     A1_SGL{i} = ImpactFunc( model_SGL(i), options );
%     LogLike_MLE_SGL(i) = Loglike_Basis( Seqs1(2001:3000), model_SGL(i), alg1 );
% 
% end
% 
% % Visualize the infectivity matrix (the adjacent matrix of Granger causality graph)
% for k = 1:nn
%     figure
%     subplot(221)        
%     imagesc(A)
%     title('Ground truth of infectivity')
%     axis square
%     colorbar
% 
%     subplot(222)        
%     imagesc(A1{k})
%     title(sprintf('Estimated MLE TNum=%d ', 500 + (k-1)*500))
%     colorbar
%     axis square
% 
%     subplot(223)   
%     imagesc(A1_S{k})
%     title(sprintf('Estimated MLE-S TNum=%d ', 500 + (k-1)*500))
%     colorbar
%     axis square
% 
%     subplot(224)        
%     imagesc(A1_SGL{k})
%     title(sprintf('Estimated MLE-SGL TNum=%d ', 500 + (k-1)*500))
%     colorbar
%     axis square
% 
%     NUM = 500:500:2000;
% 
% 
% end
%     
%     figure 
%     hold on
%     plot(NUM, LogLike_MLE_SGL(:), 'p-');
%     plot(NUM, LogLike_MLE_S(:), 'g-');
%     plot(NUM, LogLike_MLE(:), 'b-'); 
%     legend( 'MLE-SGL','MLE-S','MLE');
% 
%     axis tight
%     axis square
%     ylabel('LogLike');
%     hold off
% 
