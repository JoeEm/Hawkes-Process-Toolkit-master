%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Learning Granger causality for Hawkes processes by MLE based method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear


%for Tmax 1 = 1honr, 1 week = 24*7 = 168, 
switchFlag = 1;

options.N = 500; % the number of sequences                       2000
options.Nmax = 500; % the maximum number of events per sequence 500
options.Tmax = 500; % the maximum size of time window            200
options.tstep = 0.1;
options.dt = 0.1; % the length of each time step
options.M = 250; % the number of steps in the time interval for computing sup-intensity
options.GenerationNum = 100; % the number of generations for branch processing
D = 8; % the dimension of Hawkes processes


disp('Approximate simulation of Hawkes processes via branching process')
[Seqs1,options] = GeneratingSimulationData(options,D);


disp('Maximum likelihood estimation and basis representation') 

alg1.LowRank = 0; % without low-rank regularizer
alg1.Sparse = 1; % with sparse regularizer
alg1.alphaS = 10; %MLE-S
alg1.GroupSparse = 1; % with group-sparse regularizer
alg1.alphaGS = 100; %MLE-SGL
alg1.outer = 5;%5
alg1.rho = 0.1; % the initial parameter for ADMM
alg1.inner = 8;%8
alg1.thres = 1e-5;
alg1.Tmax = [];
alg1.storeErr = 0;
alg1.storeLL = 0;
alg1.As = 10; %for Local Independence R.
alg1.Ap = 1000;%for Pairwise similarity R.



%model1 = Initialization_Basis(Seqs1);
% model1.kernel = 'gauss';
% model1.w = pi*options.Nmax/(options.Tmax);  %According to the Real-Data Experiment
% model1.landmark = (model1.w)*(0:ceil(options.Tmax/model1.w));
% model1.A = rand(D, length(model1.landmark), D);
% model1.mu = rand(D, 1)./D;


%learning the model by MLE (first half )
[NewSeqs1, NewSeqs2]  = BreakSeqs(Seqs1);
model1 = Initialization_Basis(NewSeqs1); 
model1_A = Learning_MLE_Basis( NewSeqs1, model1, alg1 ); 
model1_result_A = ImpactFunc( model1_A, options );

%learning the model by MLE (second half )
model2 = Initialization_Basis(NewSeqs2); 
model2_A = Learning_MLE_Basis( NewSeqs2, model2, alg1 ); 
model2_result_A = ImpactFunc( model2_A, options );



% Visualize the infectivity matrix (the adjacent matrix of Granger causality graph)
figure
subplot(121)        
imagesc(model1_result_A)
title('Estimated infectivity-first half') 
axis square
colorbar
subplot(122)        
imagesc(model2_result_A)
title('Estimated infectivity-second-half')
colorbar
axis square
