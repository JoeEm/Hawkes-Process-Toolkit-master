%clear;

% options.N = 500; % the number of sequences
% options.Nmax = 500; % the maximum number of events per sequence
% options.Tmax = 200; % the maximum size of time window
% options.tstep = 0.1;
% options.dt = 0.1; % the length of each time step
% options.M = 250; % the number of steps in the time interval for computing sup-intensity
% options.GenerationNum = 100; % the number of generations for branch processing
% D = 8; % the dimension of Hawkes processes
% 
% Group = 2;
% for k = 1:Group
%     [Seqs{k},para] = GenerateSingleGroupData(options,D,k);
%     % Visualize all impact functions and infectivity matrix
%     %first test: only exsit two patterns.
%     A{k} = ImpactFunc( para, options );
% end
% 
% figure
% subplot(121)        
% imagesc(A{1})
% title('Ground truth of infectivity G1')
% axis square
% colorbar
% 
% subplot(122)        
% imagesc(A{2})
% title('Ground truth of infectivity G2')
% axis square
% colorbar
% 
% Seqs = LinkingTwoSeqsToOne(Seqs{1},Seqs{2});
% 
% 
% Seqs1 = Seqs;
% 
% disp('Maximum likelihood estimation and basis representation') 
% 
% NewSeqs  = BreakSeqs(Seqs1,options);
% NewSeqs1 = NewSeqs{1};
% NewSeqs2 = NewSeqs{2};


alg1.LowRank = 0; % without low-rank regularizer
alg1.Sparse = 1; % with sparse regularizer
alg1.alphaS = 10; %MLE-S
alg1.GroupSparse = 1; % with group-sparse regularizer
alg1.alphaGS = 100; %MLE-SGL
alg1.outer = 5;%5
alg1.rho = 0.1; % the initial parameter for ADMM--
alg1.inner = 8;%8
alg1.thres = 1e-5;
alg1.Tmax = [];
alg1.storeErr = 0;
alg1.storeLL = 0;
alg1.As = 10; %for Local Independence R.
alg1.Ap = 1000;%for Pairwise similarity R.


model1 = Initialization_Basis(Seqs1{1}); 
model2 = Initialization_Basis(Seqs1{2}); 

%model1 = Initialization_Basis(Seqs1);
% model1.kernel = 'gauss';
% model1.w = pi*options.Nmax/(options.Tmax);  %According to the Real-Data Experiment
% model1.landmark = (model1.w)*(0:ceil(options.Tmax/model1.w));
% model1.A = rand(D, length(model1.landmark), D);
% model1.mu = rand(D, 1)./D;


% learning the model by MLE-SGL
% alg1.Sparse=1;
% alg1.GroupSparse=1; 
% model_MLE_SGL_1 = Learning_MLE_S_Basis( Seqs1{1}, model1, alg1 );
% model_MLE_SGL_2 = Learning_MLE_S_Basis( Seqs1{2}, model2, alg1 );
% ClusterNumbers = 5;
% for i = 1:ClusterNumbers
%     
%     para.kernel = paraDataGenrating{i}.kernel;
%     para.w = paraDataGenrating{i}.w;
%     para.landmark = paraDataGenrating{i}.landmark;
%     para.mu = paraDataGenrating{i}.mu;
%     para.A = paraDataGenrating{i}.A;
%     options.Tmax = 200;
%     A{i} = ImpactFunc(para, options);
% end
a = 9;

ClusterNumbers = 5;
for i = 1:ClusterNumbers
    para.kernel = modelStorage{a}(i).kernel;
    para.w = modelStorage{a}(i).w;
    para.landmark = modelStorage{a}(i).landmark;
    para.mu = modelStorage{a}(i).mu;
    para.A = modelStorage{a}(i).A;
    options.Tmax = 200;
    A{i} = ImpactFunc(para, options);
end


figure

subplot(231)        
imagesc(A{1})
%title(sprintf('Ground-Truth-MLE First-Half'))
title(sprintf('Estimated-MLE-1')) 
axis square
colorbar

subplot(232)       
imagesc(A{2})
%title(sprintf('Ground-Truth-MLE Second-Half'))
title(sprintf('Estimated-MLE-2'))  
colorbar
axis square

subplot(233)        
imagesc(A{3})
%title(sprintf('Ground-Truth-MLE First-Half'))
title(sprintf('Estimated-MLE-3'))   
axis square
colorbar

subplot(234)        
imagesc(A{4})
%title(sprintf('Ground-Truth-MLE First-Half')) 
title(sprintf('Estimated-MLE-4'))  
axis square
colorbar

subplot(235)        
imagesc(A{5})
%title(sprintf('Ground-Truth-MLE First-Half')) 
title(sprintf('Estimated-MLE-5'))  
axis square
colorbar

% 
% A1 = ImpactFunc(model_MLE_SGL, options);
% 
% 
% model_MLE_SGL = Learning_MLE_S_Basis( NewSeqs2, model2, alg1 );
% 
% A2 = ImpactFunc(model_MLE_SGL, options);


% Visualize the infectivity matrix (the adjacent matrix of Granger causality graph)

% for i = 1:5 
%     for j = 1:2
%         options.kernel = modelStorage{i}(j).kernel;
%         options.w = modelStorage{i}(j).w;
%         options.landmark = modelStorage{i}(j).landmark;
%         options.mu = modelStorage{i}(j).mu;
%         
%         A{i}{j} = ImpactFunc(modelStorage{i}(j), options);
%     end
% end
% 
% for i = 1:5 
%     figure
%     subplot(121)        
%     imagesc(A{i}{1})
%     title(sprintf('Estimated infectivity-MLE First-Half_{%d}',i)) 
%     axis square
%     colorbar
%     subplot(122)       
%     imagesc(A{i}{2})
%     title(sprintf('Estimated infectivity-MLE Second-Half_{%d}',i)) 
%     colorbar
%     axis square
% 
% end

