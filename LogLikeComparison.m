%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Learning Granger causality for Hawkes processes by MLE based method
%
% IPTVData Exeperiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear

 
options.N = 2000; % the number of sequences   2000
options.Nmax = 500; % the maximum number of events per sequence when the N = 200 (knda Sparse)
options.Tmax = 24*8; % according to the Real-Data Experiment (time length is 8 days)
options.tstep = 0.1;
options.dt = 0.1; % the length of each time step
options.M = 250; % the number of steps in the time interval for computing sup-intensity
options.GenerationNum = 100; % the number of generations for branch processing
D = 16; % the dimension of Hawkes processes



disp('Approximate simulation of Hawkes processes via branching process')
disp('Complicated gaussian kernel')
para1.kernel = 'gauss'; % the type of kernels per impact function
para1.w = 2; % the bandwidth of gaussian kernel
para1.landmark = 0:4:12; % the central locations of kernels
L = length(para1.landmark);

% initialize ground truth parameters
% para1.mu = rand(D,1)/D;
% para1.A = zeros(D, D, L);
% for l = 1:L
%     para1.A(:,:,l) = (0.5^l)*(0.5+ones(D));
% end
% mask = rand(D).*double(rand(D)>0.7);
% para1.A = para1.A.*repmat(mask, [1,1,L]);
% para1.A = 0.25*para1.A./max(abs(eig(sum(para1.A,3)))); % ensure the stationarity of Hawkes process
% tmp = para1.A;
% para1.A = reshape(para1.A, [D, L, D]);
% 
% for di = 1:D
%     for dj = 1:D
%         phi = tmp(di, dj, :);
%         para1.A(dj, :, di) = phi(:);
%     end
% end
% 
% for di = 1:D
%     for dj = 1:D
%         if (di == 1|| di == 2 || di == 3 || dj ==1 || dj==2 || dj ==3 )
%             para1.A(dj, :, di) = 0;
%         end
%     end
% end

%do operation on A 


% two simulation methods
% Seqs1 = Simulation_Branch_HP(para1, options);
% Seqs1 = Simulation_Thinning_HP(para1, options);

load IPTVData.mat
%Seqs1 = Seqs(1:300);
Seqs1 = Seqs;
Seqs1 = Seqs(1:30);
%Seqs1 = [Seqs1, Seqs(54:91)];
Seqs1 = [Seqs1, Seqs(137:173)];
%Seqs1 = [Seqs1, Seqs(249:271)];

NewSeqs = IPTV_RealData_Cutting(Seqs1);


alg1.LowRank = 0; % without low-rank regularizer
alg1.Sparse = 1; % with sparse regularizer
alg1.alphaS = 10; %MLE-S
alg1.GroupSparse = 1; % with group-sparse regularizer
alg1.alphaGS = 100; %MLE-SGL
alg1.outer = 5;
alg1.rho = 0.1; % the initial parameter for ADMM
alg1.inner = 8;
alg1.thres = 1e-5;
alg1.Tmax = [];
alg1.storeErr = 0;
alg1.storeLL = 0;
alg1.As = 10; %for Local Independence R.
alg1.Ap = 1000;%for Pairwise similarity R.





nn = 4; %5iteration
for i = 1:nn
    %num = 250 + 50*i;
    %num = 30;
    
    %Seqs1 = ExtractingIPTVData(Seqs1, i);
    OutPutSeqs = RearrangeTheSeqs(NewSeqs, i)
    TraingingSampleNum = length(OutPutSeqs);
    Seqs1 =  OutPutSeqs;
    
    sample = 0;
    for k = 1:TraingingSampleNum
        sample = sample + length(OutPutSeqs(k).Time);
    end
    sample = sample/TraingingSampleNum;
    
    if i == 1
       %model1 = Initialization_Basis(Seqs1);
       model1.kernel = 'gauss';
       model1.w = pi*sample/(24*7);  %According to the Real-Data Experiment
       model1.landmark = (model1.w)*(0:ceil(options.Tmax/model1.w));
       model1.A = rand(D, length(model1.landmark), D);
       model1.mu = rand(D, 1)./D;

       ImpOptions.landmark = model1.landmark;
       ImpOptions.Tmax = options.Tmax;
       A_GT = ImpactFunc( model1, ImpOptions);
    end
    disp('learning the model by MLE'); 
    % learning the model by MLE
    alg1.Sparse=0;
    alg1.GroupSparse=0; 
    %alg1.LowRank=0;
    %model_MLE(i) = Learning_MLE_Basis( Seqs1(1:(50*i)), model1, alg1 );
    %model_MLE(i) = Learning_MLE_Basis( Seqs1(5*(4*(1 + 3*(i - 1))):TraingingSampleNum), model1, alg1 );
    model_MLE(i) = Learning_MLE_Basis( Seqs1(101:TraingingSampleNum), model1, alg1 );
    
    disp('calculating the RelErrMu and RelErrImFunc MLE'); 
    [RelErrMu_MLE(i), RelErrImFunc_MLE(i)] = ReErr_MuAndImpFunc( options.Tmax, model1, model_MLE(i));
    disp('Loglike of MLE'); 
    Loglike_MLE(i)= Loglike_Basis( Seqs1(1:100), model_MLE(i), alg1 );
    A_MLE{i} = ImpactFunc( model_MLE(i), ImpOptions );
    
    
    if i == 1
       %model2 = Initialization_Basis(Seqs1(6:10));
       model2 = model1;
    end
    disp('learning the model by MLE-S'); 
    % learning the model by MLE-S
    alg1.Sparse=1;
    alg1.GroupSparse=0; 
    %alg1.LowRank=0;
    %model_MLE_S(i) = Learning_MLE_S_Basis( Seqs1(1:(50*i)), model2, alg1 );
    %model_MLE_S(i) = Learning_MLE_S_Basis(  Seqs1(5*(4*(1 + 3*(i - 1))):TraingingSampleNum) , model2, alg1 );
    model_MLE_S(i) = Learning_MLE_S_Basis( Seqs1(101:TraingingSampleNum), model1, alg1 );
  
    disp('calculating the RelErrMu and RelErrImFunc MLE-S'); 
    [RelErrMu_MLE_S(i), RelErrImFunc_MLE_S(i)] = ReErr_MuAndImpFunc( options.Tmax, model2, model_MLE_S(i));
    disp('Loglike of MLE-S'); 
    Loglike_MLE_S(i)= Loglike_Basis(Seqs1(1:100), model_MLE_S(i), alg1 );
    A_MLE_S{i} = ImpactFunc( model_MLE_S(i), ImpOptions );
    
    if i == 1
        %model3 = Initialization_Basis(Seqs1(11:15));
        model3 = model1;
    end
    disp('learning the model by MLE-SGL'); 
    % learning the model by MLE-SGL
    alg1.Sparse=1;
    alg1.GroupSparse=1; 
    %alg1.LowRank=1;
    %model_MLE_SGL(i) = Learning_MLE_S_Basis( Seqs1(1:(50*i)), model3, alg1 );
    %model_MLE_SGL(i) = Learning_MLE_S_Basis( Seqs1(5*(4*(1 + 3*(i - 1))):TraingingSampleNum), model3, alg1 );
    model_MLE_SGL(i) = Learning_MLE_S_Basis( Seqs1(101:TraingingSampleNum), model1, alg1 );
  
    disp('calculating the RelErrMu and RelErrImFunc MLE-SGL'); 
    [RelErrMu_MLE_SGL(i), RelErrImFunc_MLE_SGL(i)] = ReErr_MuAndImpFunc( options.Tmax, model3, model_MLE_SGL(i));
    disp('Loglike of MLE-SGL');
    Loglike_MLE_SGL(i)= Loglike_Basis( Seqs1(1:100), model_MLE_SGL(i), alg1 );
    A_MLE_SGL{i} = ImpactFunc( model_MLE_SGL(i), ImpOptions );

%     if i == 1
%         %model4 = Initialization_Basis(Seqs1(16:20));
%         model4 = model1;
%     end
%     disp('learning the model by MLE-SGLP'); 
%     %learning the model by MLE-SGLP
%     alg1.Sparse=1;
%     alg1.GroupSparse=1; 
%     %alg1.LowRank=1;
%     model_MLE_SGLP(i) = Learning_MLE_Basis( Seqs1(1:(50*i)), model4, alg1 );
%     disp('calculating the RelErrMu and RelErrImFunc MLE-SGLP'); 
%     [RelErrMu_MLE_SGLP(i), RelErrImFunc_MLE_SGLP(i)] = ReErr_MuAndImpFunc( options.Tmax, model4, model_MLE_SGLP(i));
%     disp('Loglike of MLE-SGLP');
%     Loglike_MLE_SGLP(i)= Loglike_Basis( Seqs1(251:num), model_MLE_SGLP(i), alg1 );
%     A_MLE_SGLP{i} = ImpactFunc( model_MLE_SGLP(i), ImpOptions );

    
%     if i == 1
%        %model5 = Initialization_Basis(Seqs1(21:25));
%         model5 = model1;
%     end
%     disp('learning the model by MLE-ODE'); 
%     alg1.alpha = alg1.alphaS;
%  
%     %learning the model by MLE-ODE
%     model_MLE_ODE(i) = Learning_MLE_ODE( Seqs1(1:(50*i)), model5, alg1 );
%     disp('calculating the RelErrMu and RelErrImFunc MLE-ODE'); 
%     [RelErrMu_MLE_ODE(i), RelErrImFunc_MLE_ODE(i)] = ReErr_MuAndImpFunc( options.Tmax, model5, model_MLE_ODE(i));
%     disp('Loglike of MLE-ODE');
%     Loglike_MLE_ODE(i)= Loglike_Basis( Seqs1(251:num), model_MLE_ODE(i), alg1 );
%     A_MLE_ODE{i} = ImpactFunc( model_MLE_ODE(i), ImpOptions );
%     
     
%     if i == 1
%         %model6 = Initialization_Basis(Seqs1(26:35));
%         model6 = model1;
%     end
%     disp('learning the model by LS'); 
%     %learning the model by MLE-LS
%     model_LS(i) = Learning_LS_Discrete( Seqs1(1:(50*i)), model6, alg1 );
%     disp('calculating the RelErrMu and RelErrImFunc LS'); 
%     [RelErrMu_LS(i), RelErrImFunc_LS(i)] = ReErr_MuAndImpFunc( options.Tmax, model1, model_LS(i));
%     disp('Loglike of LS');
%     Loglike_LS(i)= Loglike_Basis( Seqs1(251:num), model_LS(i), alg1 );
%     A_LS{i} = ImpactFunc( model_LS(i), ImpOptions );
  
end

NUM = 1:3:10;
figure   
subplot(131); %ploting the Loglike
hold on
%plot(NUM, Loglike_MLE_SGLP(:), 'r-');
plot(NUM, Loglike_MLE_SGL(:), 'p-');
plot(NUM, Loglike_MLE_S(:), 'g-');
plot(NUM, Loglike_MLE(:), 'b-');
% plot(NUM, Loglike_MLE_ODE(:), 'k-');
% plot(NUM, Loglike_LS(:), 'c-');
%legend('MLE-SGLP', 'MLE-SGL','MLE-S','sMLE','ODE','LS');
axis tight
axis square
ylabel('LogLike');
xlabel('The number of training samples');
hold off

% subplot(132); %ploting the RelErrMu
% hold on
% %plot(NUM, RelErrMu_MLE_SGLP(:), 'r-');
% plot(NUM, RelErrMu_MLE_SGL(:), 'p-');
% plot(NUM, RelErrMu_MLE_S(:), 'g-');
% plot(NUM, RelErrMu_MLE(:), 'b-');
% plot(NUM, RelErrMu_MLE_ODE(:), 'k-');
% plot(NUM, RelErrMu_LS(:), 'c-');
% %legend('MLE-SGLP', 'MLE-SGL','MLE-S','MLE','ODE','LS');
% axis tight
% axis square
% ylabel('Relative error of mu');
% xlabel('The number of training samples');
% hold off
% 
% subplot(133); %ploting the RelErrMu
% hold on
% %plot(NUM, RelErrImFunc_MLE_SGLP(:), 'r-');
% plot(NUM, RelErrImFunc_MLE_SGL(:), 'p-');
% plot(NUM, RelErrImFunc_MLE_S(:), 'g-');
% plot(NUM, RelErrImFunc_MLE(:), 'b-');
% plot(NUM, RelErrImFunc_MLE_ODE(:), 'k-');
% plot(NUM, RelErrImFunc_LS(:), 'c-');
legend( 'MLE-SGL','MLE-S','MLE');
% axis tight
% axis square
% ylabel('Relative error of A');
% xlabel('The number of training samples');
% hold off

    j = 1;

for i = 1:4
    figure

%     
%     subplot(332)        
%     imagesc(A_GT)
%     title(sprintf('Ground Truth TrainNum = %d months', 1 + (i-1)*4));
%     axis square
%     colorbar
    
    subplot(334)        
    imagesc(A_MLE{j})
    title(sprintf('Estimated MLE TrainNum = %d months', 1 + (i-1)*3 ))
    axis square
    colorbar
    
    
    subplot(335)        
    imagesc(A_MLE_S{j})
    title(sprintf('Estimated MLE-S TrainNum = %d months', 1 + (i-1)*3))
    colorbar
    axis square
    
    subplot(336)        
    imagesc(A_MLE_SGL{j})
    title(sprintf('Estimated MLE-SGL TrainNum = %d months', 1 + (i-1)*3))
    colorbar
    axis square
    
%     subplot(337)        
%     imagesc(A_MLE_SGLP{j})
%     title(sprintf('Estimated MLE-SGLP TrainNum = %d', 50 + (i - 1)*50))
%     colorbar
%     axis square
    
%     subplot(338)        
%     imagesc(A_MLE_ODE{j})
%     title(sprintf('Estimated MLE-ODE TrainNum = %d', 50 + (i - 1)*50))
%     colorbar
%     axis square
    
%     subplot(339)        
%     imagesc(A_LS{j})
%     title(sprintf('Estimated LS TrainNum = %d', 50 + (i - 1)*50))
%     colorbar
%     axis square
    
    j = j + 1;
    
end